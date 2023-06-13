use crate::genome::GenomeInfo;
use crate::vcf::*;

#[test]
fn compare_table_matrix_enc() {
    use bstr::ByteSlice;

    let genome_info = "testdata/genome.toml";
    let ginfo = GenomeInfo::from_toml_file(genome_info);
    let vcf_path = "testdata/bcf/small.bcf";

    let max_maf = 0.001f64;
    let (sit, _ind, rec) = read_vcf(&ginfo, vcf_path, max_maf, None);

    // here use a very low min_maf to encode rare variants and common varints in matrix,
    // so that it can use be to very the accuracy of table encoding of rare variants
    let min_maf = 0.0;
    let (sit2, _ind2, gm2) = read_vcf_for_genotype_matrix(&ginfo, vcf_path, min_maf, None);

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
