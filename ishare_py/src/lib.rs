use ahash::AHashSet;
use ahash::HashMap;
use ahash::HashMapExt;

use ishare::vcf::read_vcf;
use ishare::{
    genome::GenomeInfo,
    genotype::{afreq::get_afreq_from_vcf_genome_wide, rare::GenotypeRecords},
    gmap::GeneticMap,
    indiv::Individuals,
    share::ibd::{
        coverage::CovCounter,
        ibdseg::IbdSeg,
        ibdset::{IbdSet, IbdSetPloidyStatus, IbdSetSortStatus},
    },
    site::Sites,
    stat::xirs::{XirsBuilder, XirsBuilder2},
};
use numpy::IntoPyArray;
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyDict, PyList},
};

#[pyclass]
struct PyGeneticMap {
    gmap: GeneticMap,
}
#[pymethods]
impl PyGeneticMap {
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
            let gmap: &PyCell<PyGeneticMap> = gmap.extract()?;
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
        Ok(PyGeneticMap {
            gmap: GeneticMap::from_vec(out),
        })
    }
}

/// a version of IbdSet to construct pyclass (no references)
#[pyclass]
struct PyIbdSet {
    gmap: GeneticMap,
    ginfo: GenomeInfo,
    inds: Individuals,
    ibd: Vec<IbdSeg>,
    ploidy_status: IbdSetPloidyStatus,
    sort_status: IbdSetSortStatus,
}

#[pymethods]
impl PyIbdSet {
    #[new]
    fn new(
        chromsizes: Vec<u32>,
        chromnames: Vec<String>,
        samples: Vec<String>,
        gmaps: Vec<Vec<(u32, f32)>>,
        ibd_files: Vec<String>,
        ibd_format: String,
    ) -> PyResult<Self> {
        let errm = |m| Err(PyErr::new::<PyValueError, _>(m));
        if chromsizes.len() != chromnames.len() {
            return errm("chromsizes length not equal to chrnames length");
        }
        if chromsizes.len() != gmaps.len() {
            return errm("chromsizes length not equal to gmaps length");
        }
        if !((chromsizes.len() == ibd_files.len()) || (ibd_files.len() == 1)) {
            return errm("ibd_files length is neither 1 nor chromsizes length");
        }
        // genetic maps
        let gmaps = GeneticMap::from_gmap_vec(&gmaps, &chromsizes);

        // genome info
        let name = "".to_owned();
        let idx: HashMap<String, usize> = chromnames
            .iter()
            .enumerate()
            .map(|(i, n)| (n.to_owned(), i))
            .collect();
        let mut r: u32 = 0;
        let mut gwstarts = Vec::<u32>::with_capacity(chromsizes.len());
        for l in chromsizes.iter() {
            gwstarts.push(r);
            r += l;
        }
        let ginfo = GenomeInfo::new_from_parts(name, chromsizes, chromnames.clone(), idx, gwstarts);

        // individuals
        let inds = Individuals::from_iter(samples.iter().map(|s| s.as_str()));

        // Ibdset
        let mut ibdset = IbdSet::new(&gmaps, &ginfo, &inds);

        for (fname, chrname) in ibd_files.iter().zip(chromnames.iter()) {
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

        Ok(PyIbdSet {
            gmap: gmaps,
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
struct RareVariants {
    records: GenotypeRecords,
    sort_status: u8,
    inds: Option<Individuals>,
    sits: Option<Sites>,
}

#[pymethods]
impl RareVariants {
    #[staticmethod]
    fn from_vcf<'py>(
        chrnames: Vec<String>,
        chrsizes: Vec<u32>,
        samples: Vec<String>,
        vcf: &str,
        chunksize: usize,
        max_maf: f64,
        _py: Python<'py>,
    ) -> PyResult<Self> {
        // create ginfo
        let mut gwstarts = vec![0];
        let mut idx = HashMap::<String, usize>::new();
        for (i, chrlen) in chrsizes.iter().enumerate() {
            idx.insert(chrnames[i].to_owned(), i);
            if gwstarts.len() < chrsizes.len() {
                let last = gwstarts.last().unwrap();
                gwstarts.push(last + chrlen);
            }
        }
        let ginfo =
            GenomeInfo::new_from_parts("genome".to_owned(), chrsizes, chrnames, idx, gwstarts);

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

        Ok(RareVariants {
            records,
            sort_status: 1u8,
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

        self.records.clone().into_parquet_file(&gt_file);
        self.sits.clone().unwrap().into_parquet_file(&sites_file);
        self.inds.clone().unwrap().into_parquet_file(&ind_file);
        Ok(())
    }

    #[staticmethod]
    fn from_files(prefix: String, read_sites: bool, read_individuals: bool) -> PyResult<Self> {
        use ishare::utils::path::from_prefix;
        let gt_file = from_prefix(&prefix, "rec").unwrap();
        let sites_file = from_prefix(&prefix, "sit").unwrap();
        let ind_file = from_prefix(&prefix, "ind").unwrap();
        let records = GenotypeRecords::from_parquet_file(&gt_file);
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
        Ok(RareVariants {
            records,
            sort_status: 0, // TODO: need check
            inds,
            sits,
        })
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn isharepy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIbdSet>()?;
    m.add_class::<PyGeneticMap>()?;
    Ok(())
}
