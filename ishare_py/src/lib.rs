use ahash::{AHashSet, HashMap, HashMapExt};

use ishare::{
    genome::GenomeInfo,
    genotype::{afreq::get_afreq_from_vcf_genome_wide, rare::GenotypeRecords},
    gmap::GeneticMap,
    indiv::Individuals,
    share::{
        ibd::{
            coverage::CovCounter,
            ibdseg::IbdSeg,
            ibdset::{IbdSet, IbdSetPloidyStatus, IbdSetSortStatus},
        },
        mat::NamedMatrix,
        pairs::PairChunkIter,
    },
    site::Sites,
    stat::xirs::{XirsBuilder, XirsBuilder2},
    utils::path::from_prefix,
    vcf::read_vcf,
};
use numpy::{IntoPyArray, ToPyArray};
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyDict, PyList},
};
use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

#[pyclass]
struct GInfo {
    ginfo: GenomeInfo,
}

#[pymethods]
impl GInfo {
    #[new]
    fn new(
        name: String,
        chromsizes: Vec<u32>,
        chromnames: Vec<String>,
        additional_chromname_2_id_map: HashMap<String, usize>,
    ) -> PyResult<Self> {
        let errm = |m| Err(PyErr::new::<PyValueError, _>(m));
        if chromsizes.len() != chromnames.len() {
            return errm("chromsizes length not equal to chrnames length");
        }
        // add chromnames and its order into idx map
        let mut idx: HashMap<String, usize> = chromnames
            .iter()
            .enumerate()
            .map(|(i, n)| (n.to_owned(), i))
            .collect();
        // validate the length of chromnames
        if idx.len() != chromsizes.len() {
            return errm("the chromnames should be unique and has the same length as chromsizes");
        }
        // add additional name to id mappings
        for (k, v) in additional_chromname_2_id_map.iter() {
            match idx.get(k) {
                Some(v2) => {
                    if v != v2 {
                        return errm("chromname conflicts between chromnames and additional_chromname_2_id_map");
                    }
                }
                None => {
                    idx.insert(k.to_owned(), *v);
                }
            }
        }

        // calculate cumsum for gw_starts
        let gwstarts = {
            let mut s = 0;
            let mut v = vec![0u32];
            for x in chromsizes.iter().take(chromsizes.len() - 1) {
                s += *x;
                v.push(s);
            }
            v
        };
        let ginfo = GenomeInfo::new_from_parts(name, chromsizes, chromnames, idx, gwstarts);

        Ok(GInfo { ginfo })
    }
}

#[pyclass]
struct GMap {
    gmap: GeneticMap,
}
#[pymethods]
impl GMap {
    #[staticmethod]
    /// position is 0-based, position indicate the left bound of a base-pair
    fn from_list(mut lst: Vec<(u32, f32)>, chrlen: u32) -> PyResult<Self> {
        // trimming
        lst.retain(|x| x.0 < chrlen);
        if lst.len() < 2 {
            return Err(PyErr::new::<PyValueError, _>(
                "argument `lst`  must have a length >=2 after trimming ends",
            ));
        }

        // add left ends
        if lst.first().unwrap().0 != 0 {
            lst.insert(0, (0, 0.0));
        }

        // add right end
        let (last_bp, last_cm) = lst.last().unwrap();
        let (last_bp, last_cm) = (*last_bp, *last_cm);
        let average_rate = last_cm / (last_bp + 1) as f32;
        if last_bp != chrlen - 1 {
            lst.push((chrlen - 1, (chrlen - 1) as f32 * average_rate));
        }

        let gmap = GeneticMap::from_vec(lst);
        Ok(Self { gmap })
    }
    #[staticmethod]
    fn from_plink_map_file(p: &str, chrlen: u32) -> PyResult<Self> {
        Ok(Self {
            gmap: GeneticMap::from_plink_map(p, chrlen),
        })
    }
    #[staticmethod]
    fn from_constant_rate(rate: f32, chrlen: u32) -> PyResult<Self> {
        Ok(Self {
            gmap: GeneticMap::from_constant_rate(rate, chrlen),
        })
    }

    #[staticmethod]
    fn from_list_of_gmap<'py>(lst_gmap: &'py PyList, _py: Python<'py>) -> PyResult<Self> {
        let (mut bp_offset, mut cm_offset) = (0, 0.0);
        let mut out = vec![];
        for gmap in lst_gmap.into_iter() {
            let gmap: &PyCell<GMap> = gmap.extract()?;
            let gmap = gmap.borrow();
            let v = gmap.gmap.as_slice();
            // add data
            for (chr_bp, chr_cm) in v {
                out.push((chr_bp + bp_offset, chr_cm + cm_offset))
            }

            // update offsets
            let (bp1, cm1) = v[v.len() - 1];
            let (bp2, cm2) = v[v.len() - 2];
            let last_rate = (cm1 - cm2) / (bp1 - bp2) as f32;
            let bp_len = bp1 + 1;
            let cm_len = cm1 + last_rate * 1.0;
            bp_offset += bp_len;
            cm_offset += cm_len;
        }
        Ok(GMap {
            gmap: GeneticMap::from_vec(out),
        })
    }
}

/// a version of IbdSet to construct pyclass (no references)
#[pyclass]
struct IBD {
    gmap: GeneticMap,
    ginfo: GenomeInfo,
    inds: Individuals,
    ibd: Vec<IbdSeg>,
    ploidy_status: IbdSetPloidyStatus,
    sort_status: IbdSetSortStatus,
}

#[pymethods]
impl IBD {
    #[new]
    fn new(
        ginfo: &GInfo,
        gmap: &GMap,
        ibd_files: Vec<String>,
        ibd_format: String,
        samples: Vec<String>,
    ) -> PyResult<Self> {
        // individuals
        let inds = Individuals::from_iter(samples.iter().map(|s| s.as_str()));
        let ginfo = ginfo.ginfo.clone();
        let gmap = gmap.gmap.clone();

        // Ibdset
        let mut ibdset = IbdSet::new(&gmap, &ginfo, &inds);

        for (fname, chrname) in ibd_files.iter().zip(ginfo.chromnames.iter()) {
            if ibd_format == "hapibd" {
                ibdset.read_hapibd_file(fname, None);
            } else if ibd_format == "tskibd" {
                ibdset.read_tskibd_file(fname, chrname);
            } else if ibd_format == "hmmibd" {
                ibdset.read_hmmibd_file(fname, None);
            } else {
                return Err(PyErr::new::<PyValueError, _>(
                    "ibd_format is not correct: should be in 'hapibd', 'tskibd', or 'hmmibd'",
                ));
            }
        }

        ibdset.sort_by_haplotypes();
        ibdset.infer_ploidy();

        let (ibd, ps, ss) = ibdset.into_parts();

        Ok(IBD {
            gmap,
            ginfo,
            inds,
            ibd,
            ploidy_status: ps,
            sort_status: ss,
        })
    }

    fn get_coverage<'py>(&self, py: Python<'py>, step: u32) -> PyResult<&'py PyDict> {
        let genomesize = self.ginfo.get_total_len_bp();
        let mut cov_counter = CovCounter::from_range(0, genomesize, step);
        for seg in self.ibd.iter() {
            cov_counter.count_over_interval(&(seg.s..seg.e));
        }

        let n_sample_points = cov_counter.get_n_interval();
        let mut gw_bp = Vec::<u32>::with_capacity(n_sample_points);
        let mut ch_bp = Vec::<u32>::with_capacity(n_sample_points);
        let mut chrid = Vec::<u32>::with_capacity(n_sample_points);
        let mut cov = Vec::<u32>::with_capacity(n_sample_points);
        for (bp, _, cnt) in cov_counter.iter_sorted_start_end_count() {
            gw_bp.push(bp);
            let (chid0, _, chpos0) = self.ginfo.to_chr_pos(bp);
            chrid.push(chid0 as u32);
            ch_bp.push(chpos0);
            cov.push(cnt as u32);
        }

        let dict = PyDict::new(py);
        dict.set_item("ChromosomeId", chrid.into_pyarray(py))?;
        dict.set_item("ChrPos", ch_bp.into_pyarray(py))?;
        dict.set_item("GwPos", gw_bp.into_pyarray(py))?;
        dict.set_item("Coverage", cov.into_pyarray(py))?;
        Ok(dict)
    }

    fn get_xirs<'py>(
        &mut self,
        py: Python<'py>,
        vcf_files: Vec<String>,
        min_maf: f32,
        method: u32,
    ) -> PyResult<&'py PyDict> {
        // get allele frequency and site postion vectors from the vec files
        let mut pos_afreq_vec = get_afreq_from_vcf_genome_wide(&vcf_files, &self.ginfo);

        // remove sites of low minor allele frequency
        pos_afreq_vec.retain(|x| (x.1 >= min_maf) && (x.1 <= 1.0 - min_maf));

        let afrq: Vec<_> = pos_afreq_vec
            .iter()
            .map(|(_pos, afreq)| *afreq as f64)
            .collect();
        let site_pos: Vec<_> = pos_afreq_vec.iter().map(|(pos, _afreq)| *pos).collect();

        let mut ibdvec = vec![];
        std::mem::swap(&mut self.ibd, &mut ibdvec);
        let mut ibdset = IbdSet::new(&self.gmap, &self.ginfo, &self.inds);
        ibdset.update_parts(ibdvec, self.ploidy_status, self.sort_status);
        dbg!(self.ploidy_status);

        // calculate xirs stats
        let xirs_res = match method {
            1 => XirsBuilder::new(afrq, site_pos, &ibdset).finish(),
            2 => XirsBuilder2::new(afrq, site_pos, &ibdset).finish(),
            _ => {
                return Err(PyErr::new::<PyValueError, _>(format!(
                    "method: {method} is not implemented"
                )));
            }
        };

        // return ibd vec back to PyIbdSet
        self.ibd = ibdset.into_parts().0;

        let res_dict: &PyDict = PyDict::new(py);
        res_dict.set_item("ChrId", xirs_res.chr_id.into_pyarray(py))?;
        res_dict.set_item("ChrPos", xirs_res.chr_pos.into_pyarray(py))?;
        res_dict.set_item("GwPos", xirs_res.gw_pos.into_pyarray(py))?;
        res_dict.set_item("Raw", xirs_res.raw.into_pyarray(py))?;
        res_dict.set_item("Xirs", xirs_res.xirs.into_pyarray(py))?;
        res_dict.set_item("Pvalue", xirs_res.pval.into_pyarray(py))?;

        Ok(res_dict)
    }

    fn cmp_with(&self) {}
}

#[pyclass]
struct RVar {
    records: Option<GenotypeRecords>,
    inds: Option<Individuals>,
    sits: Option<Sites>,
}

#[pymethods]
impl RVar {
    #[staticmethod]
    fn from_vcf<'py>(
        ginfo: &GInfo,
        samples: Vec<String>,
        vcf: &str,
        chunksize: usize,
        max_maf: f64,
        _py: Python<'py>,
    ) -> PyResult<Self> {
        let ginfo = &ginfo.ginfo;

        // split chromosome into regions for parallele reading
        // divide genome into 10Mb chunks
        let regions = ginfo.partition_genome(Some(chunksize as u32));

        let target_samples: AHashSet<String> = samples.into_iter().collect();

        use rayon::prelude::*;
        // parallel running
        let mut res: Vec<_> = regions
            .into_par_iter()
            .map(|region| {
                // (sites, individuals, GenotypeRecords)
                let res = read_vcf(&target_samples, &ginfo, vcf, max_maf, region);
                if region.is_some() {
                    println!("{:?}", region);
                }
                res
            })
            .collect();

        // merge results
        let (mut sites, individuals, mut records) = res.pop().unwrap();

        for (ss, _, rr) in res {
            sites.merge(ss);
            records.merge(rr);
        }

        // sort
        _ = sites.sort_by_position_then_allele();
        records.sort_by_genome();

        Ok(RVar {
            records: Some(records),
            inds: Some(individuals),
            sits: Some(sites),
        })
    }

    fn to_files(&self, prefix: String) -> PyResult<()> {
        // construct output file names
        use ishare::utils::path::from_prefix;
        let gt_file = from_prefix(&prefix, "rec").unwrap();
        let sites_file = from_prefix(&prefix, "sit").unwrap();
        let ind_file = from_prefix(&prefix, "ind").unwrap();

        match self.records.as_ref() {
            Some(records) => {
                records.clone().into_parquet_file(&gt_file);
            }
            None => {}
        }
        match self.sits.as_ref() {
            Some(sits) => {
                sits.clone().into_parquet_file(&sites_file);
            }
            None => {}
        }
        match self.inds.as_ref() {
            Some(inds) => {
                inds.clone().into_parquet_file(&ind_file);
            }
            None => {}
        }
        Ok(())
    }

    #[new]
    fn new() -> Self {
        Self {
            records: None,
            inds: None,
            sits: None,
        }
    }

    fn read_genotype(&mut self, prefix: String) -> PyResult<()> {
        let gt_file = from_prefix(&prefix, "rec")
            .map_err(|_e| PyErr::new::<PyValueError, _>("err in path to genotype file"))?;
        self.records = Some(GenotypeRecords::from_parquet_file(gt_file));
        Ok(())
    }

    fn read_sites(&mut self, prefix: String) -> PyResult<()> {
        let sites_file = from_prefix(&prefix, "sit")
            .map_err(|_e| PyErr::new::<PyValueError, _>("err in path to sites file"))?;
        self.sits = Some(Sites::from_parquet_file(sites_file));
        Ok(())
    }

    fn read_inds(&mut self, prefix: String) -> PyResult<()> {
        let inds_file = from_prefix(&prefix, "ind")
            .map_err(|_e| PyErr::new::<PyValueError, _>("err in path to sites file"))?;
        self.inds = Some(Individuals::from_parquet_file(inds_file));
        Ok(())
    }

    #[staticmethod]
    fn from_files(
        prefix: String,
        read_genotype: bool,
        read_sites: bool,
        read_individuals: bool,
    ) -> PyResult<Self> {
        let gt_file = from_prefix(&prefix, "rec").unwrap();
        let sites_file = from_prefix(&prefix, "sit").unwrap();
        let ind_file = from_prefix(&prefix, "ind").unwrap();
        let records = if read_genotype {
            Some(GenotypeRecords::from_parquet_file(&gt_file))
        } else {
            None
        };
        let sits = if read_sites {
            let sites = Sites::from_parquet_file(sites_file);
            Some(sites)
        } else {
            None
        };
        let inds = if read_individuals {
            let inds = Individuals::from_parquet_file(ind_file);
            Some(inds)
        } else {
            None
        };
        Ok(RVar {
            records,
            inds,
            sits,
        })
    }

    fn check_genotype(&self) -> PyResult<()> {
        if let Some(_recs) = self.records.as_ref() {
            Ok(())
        } else {
            return Err(PyErr::new::<PyValueError, _>(
                "genotype is empty - have you loaded it?",
            ));
        }
    }
    fn check_sites(&self) -> PyResult<()> {
        if let Some(_) = self.sits.as_ref() {
            Ok(())
        } else {
            return Err(PyErr::new::<PyValueError, _>(
                "sites is empty - have you loaded it?",
            ));
        }
    }
    fn check_individuals(&self) -> PyResult<()> {
        if let Some(_) = self.inds.as_ref() {
            Ok(())
        } else {
            return Err(PyErr::new::<PyValueError, _>(
                "individuals is empty - have you loaded it?",
            ));
        }
    }

    fn get_genotype_table<'py>(&self, py: Python<'py>) -> PyResult<&'py PyDict> {
        self.check_genotype()?;
        let d = PyDict::new(py);
        let recs = self.records.as_ref().unwrap();
        let n = recs.records().len();
        let mut gid = Vec::with_capacity(n);
        let mut gt = Vec::with_capacity(n);
        let mut pos = Vec::with_capacity(n);
        for r in recs.records() {
            let (p, g, a) = r.get();
            gid.push(g);
            gt.push(a);
            pos.push(p);
        }
        d.set_item("GenomeId", gid.into_pyarray(py))?;
        d.set_item("GenomeWidePosition", pos.into_pyarray(py))?;
        d.set_item("AlleleID", gt.into_pyarray(py))?;
        Ok(d)
    }
    fn get_individuals_lst(&self) -> PyResult<Vec<String>> {
        self.check_individuals()?;
        Ok(self.inds.as_ref().unwrap().v().clone())
    }
    fn get_sites_table<'py>(&self, py: Python<'py>) -> PyResult<&'py PyDict> {
        self.check_sites()?;
        let sites = self.sits.as_ref().unwrap();
        let d = PyDict::new(py);
        let p = sites.get_gw_pos_slice();
        let n = p.len();
        let mut alleles_v = vec![];
        for i in 0..n {
            let alleles_str = std::str::from_utf8(sites.get_alleles_by_idx(i))?;
            alleles_v.push(alleles_str.to_owned())
        }
        d.set_item("GenomeWidePosition", p.to_pyarray(py))?;
        d.set_item("Alleles", alleles_v)?;
        Ok(d)
    }

    fn calculate_jaccard_matrix<'py>(
        &self,
        ids1: Option<Vec<u32>>,
        ids2: Option<Vec<u32>>,
        chunksize: Option<u32>,
        copy_chunk: Option<bool>,
        py: Python<'py>,
    ) -> PyResult<&'py PyDict> {
        // check
        self.check_genotype()?;
        self.check_individuals()?;

        // normalized ids lists
        // eprintln!("normalized ids list");
        let (ids1, ids2) = match (ids1, ids2) {
            (Some(ids1), Some(ids2)) => (ids1, ids2),
            (Some(ids), None) | (None, Some(ids)) => (ids.clone(), ids),
            (None, None) => {
                let n = self.inds.as_ref().unwrap().v().len() * 2;
                let ids = (0..n as u32).collect::<Vec<u32>>();
                (ids.clone(), ids)
            }
        };
        let chunksize = chunksize.unwrap_or(10000);
        let copy_chunk = copy_chunk.unwrap_or(false);
        // chunks preparation
        let mut pair_chunk_iter =
            PairChunkIter::new(ids1.as_slice(), ids2.as_slice(), chunksize as usize);

        // let nchunk = pair_chunk_iter.get_n_chunks();
        // let mut ichunk = 0usize;
        // eprintln!("chunksize={} and there are {} chunks", chunksize, nchunk);

        // chunk temporary buffers
        let mut pairs = vec![];
        let mut related = vec![];
        let mut res_vec = vec![];
        let mut is_within = false;

        // result matrix to store results
        let mut res_matrix = NamedMatrix::new(ids1.clone(), ids2.clone());

        // iterate over pair chunks
        while pair_chunk_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within) {
            // make a copy of a subset of records to better take advantage cache-locality
            let mut gt = self.records.as_ref().unwrap();
            let store: Option<GenotypeRecords>;
            // use copy_chunk to control whether to copy chunks of genotypes or directly use the genotype
            // from self
            if copy_chunk {
                store = self
                    .records
                    .as_ref()
                    .unwrap()
                    .subset_by_genomes(&related)
                    .ok();
                gt = store.as_ref().unwrap();
            }
            // parallelization over pairs within each pair-chunks
            pairs
                .as_slice()
                .into_par_iter()
                .map(|(i, j)| {
                    let (genome1, genome2) = (*i, *j);
                    let (mut total, mut shared) = (0u32, 0u32);

                    if is_within && genome1 < genome2 {
                        for (_pos, a, b) in gt.iter_genome_pair_genotypes(genome1, genome2) {
                            match (a, b) {
                                (Some(a), Some(b)) if a == b => {
                                    shared += 1;
                                    total += 1
                                }
                                (None, None) => {}
                                (_, _) => total += 1,
                            }
                        }
                        (genome1, genome2, total, shared)
                    } else {
                        // skip calculation to avoid duplication
                        (0, 0, 0, 0)
                    }
                })
                .collect_into_vec(&mut res_vec);

            // save results from vector to matrix
            for (g1, g2, total, shared) in res_vec.iter() {
                if !((*g1 == 0) && (*g2 == 0)) {
                    // if *shared > 10 {
                    //     eprintln!("{} {} {} {}", *g1, *g2, *shared, *total);
                    // }
                    res_matrix.set_by_names(*g1, *g2, (*shared as f32) / (*total as f32));
                }
            }
            // ichunk += 1;
            // eprint!("\r {:.3} %", (ichunk * 100) as f32 / nchunk as f32);
        }
        let (row_ids, col_ids, data) = res_matrix.into_parts();
        let d = PyDict::new(py);
        let data = numpy::ndarray::Array2::from_shape_vec([row_ids.len(), col_ids.len()], data)
            .map_err(|x| PyErr::new::<PyValueError, _>(x.to_string()))?;
        d.set_item("RowIds", row_ids.into_pyarray(py))?;
        d.set_item("ColumnIds", col_ids.into_pyarray(py))?;
        d.set_item("Matrix", data.into_pyarray(py))?;
        Ok(d)
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn isharepy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<IBD>()?;
    m.add_class::<GMap>()?;
    m.add_class::<GInfo>()?;
    m.add_class::<RVar>()?;
    Ok(())
}
