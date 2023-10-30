use crate::genome::GenomeInfo;
use crate::gmap::GeneticMap;
use crate::io::IntoParquet;
use crate::share::ibd::ibdset::IbdSet;
use crate::stat::xirs::*;
use crate::vcf::*;
use ahash::AHashSet;

#[test]
fn compare_table_matrix_enc() {
    use bstr::ByteSlice;

    let genome_info = "testdata/dir001/genome.toml";
    let ginfo = GenomeInfo::from_toml_file(genome_info);
    let vcf_path = "testdata/dir001/bcf/sel_chr1.bcf";

    let max_maf = 0.001f64;
    let (sit, _ind, rec) = read_vcf(&AHashSet::new(), &ginfo, vcf_path, max_maf, None);

    // here use a very low min_maf to encode rare variants and common varints in matrix,
    // so that it can use be to very the accuracy of table encoding of rare variants
    let min_maf = 0.0;
    let (sit2, _ind2, gm2) =
        read_vcf_for_genotype_matrix(&AHashSet::new(), &ginfo, vcf_path, min_maf, None);

    for r in rec.records() {
        let pos = r.get_position();
        let aid = r.get_allele();
        let gid = r.get_genome();
        let alleles = sit.get_site_by_position(pos);

        print!(
            "alleles=|{}|, aid={aid}, pos={pos}, gid={gid}, ",
            alleles.as_bstr(),
        );
        let allele = alleles
            .as_bstr()
            .split(|x| *x == b' ')
            .nth(aid as usize)
            .unwrap();
        println!("allele=|{}|", allele.as_bstr());

        println!("sit2: {pos}, {}", sit2.get_site_by_position(pos).as_bstr());

        assert!(gm2.has_allele(pos, gid, allele, &sit2));
    }
}

#[test]
fn calc_xirs() {
    let ginfo = GenomeInfo::from_toml_file("testdata/dir001/genome.toml");
    let gmap = GeneticMap::from_genome_info(&ginfo);

    let vcf_fns = [
        "testdata/dir001/vcf_filt/sel_chr1.vcf.gz",
        "testdata/dir001/vcf_filt/sel_chr2.vcf.gz",
        "testdata/dir001/vcf_filt/sel_chr3.vcf.gz",
    ];
    let target_samples = AHashSet::new();
    let (mut sites, inds, mut mat) =
        read_vcf_for_genotype_matrix(&target_samples, &ginfo, vcf_fns[0], 0.01, None);
    for vcf_fn in &vcf_fns[1..] {
        let (sites2, inds2, mat2) =
            read_vcf_for_genotype_matrix(&target_samples, &ginfo, vcf_fn, 0.01, None);
        assert!(inds.v() == inds2.v());
        sites.merge(sites2);
        mat.merge(mat2);
    }

    println!("mat size: nrow={}, ncol={}", mat.nrows(), mat.ncols());
    // return;

    let mut ibd = IbdSet::new(&gmap, &ginfo, &inds);
    ibd.read_hapibd_dir("testdata/dir001/ibd_hapibd/");
    ibd.sort_by_haplotypes();
    ibd.infer_ploidy();

    let afreq = mat.get_afreq();
    let gw_pos = sites.get_gw_pos_slice().to_owned();

    let mut xirs = XirsBuilder::new(afreq, gw_pos, &ibd).finish();

    xirs.into_parquet("tmp.xirs.pq");
}
