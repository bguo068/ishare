use ahash::HashMap;

use ishare::{
    genome::GenomeInfo,
    gmap::GeneticMap,
    indiv::Individuals,
    share::ibd::{
        coverage::CovCounter,
        ibdseg::IbdSeg,
        ibdset::{IbdSet, IbdSetPloidyStatus, IbdSetSortStatus},
    },
    stat::xirs::{XirsBuilder, XirsBuilder2},
};
use numpy::IntoPyArray;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

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
        if chromsizes.len() != chromnames.len() {
            return Err(PyErr::new::<PyValueError, _>(
                "chromsizes length not equal to chrnames length",
            ));
        }
        if chromsizes.len() != gmaps.len() {
            return Err(PyErr::new::<PyValueError, _>(
                "chromsizes length not equal to gmaps length",
            ));
        }
        if !((chromsizes.len() == ibd_files.len()) || (ibd_files.len() == 1)) {
            return Err(PyErr::new::<PyValueError, _>(
                "ibd_files length is neither 1 nor chromsizes length",
            ));
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
        use ishare::genotype::afreq::get_afreq_from_vcf_genome_wide;
        let mut pos_afreq_vec = get_afreq_from_vcf_genome_wide(&vcf_files, &self.ginfo);
        // eprintln!("pos_len = {}", pos_afreq_vec.len());
        // remove sites of low minor allele frequency
        pos_afreq_vec.retain(|x| (x.1 >= min_maf) && (x.1 <= 1.0 - min_maf));

        // eprintln!("pos_len = {}", pos_afreq_vec.len());
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
        // let xirs_res = XirsBuilder::new(afrq, site_pos, &ibdset).finish();
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

/// A Python module implemented in Rust.
#[pymodule]
fn isharepy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIbdSet>()?;
    Ok(())
}
