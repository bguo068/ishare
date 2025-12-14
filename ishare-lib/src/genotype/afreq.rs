use ahash::{HashMap, HashMapExt};
use rust_htslib::{
    bcf::{self, record::GenotypeAllele, Read},
    htslib::bcf_is_snp,
};
use snafu::prelude::*;
use std::backtrace::Backtrace;

use crate::genome::GenomeInfo;

#[derive(Snafu, Debug)]
pub enum Error {
    HtslibError {
        source: rust_htslib::errors::Error,
        backtrace: Option<Backtrace>,
    },
    VcfMissingRid,
    Utf8Error {
        source: std::str::Utf8Error,
        backtrace: Option<Backtrace>,
    },
}

type Result<T> = std::result::Result<T, Error>;

/// similar to get_afreq_from_vcf but use genome-wide coordinates
///
/// return a vector of tuple with element being a tuple of gw_pos and afreq
pub fn get_afreq_from_vcf_genome_wide(
    vcf_files: &[String],
    ginfo: &GenomeInfo,
) -> Result<Vec<(u32, f32)>> {
    let res_map = get_afreq_from_vcf(vcf_files)?;
    let mut res_vec = Vec::<(u32, f32)>::new();
    println!("number chromosome: {}", res_vec.len());
    for (chrname, v) in res_map.iter() {
        let chrid = ginfo.idx[chrname];
        for (pos, afreq) in v.iter() {
            let gw_pos = ginfo.to_gw_pos(chrid, *pos);
            res_vec.push((gw_pos, *afreq));
        }
    }
    res_vec.sort_by_key(|tuple| tuple.0);
    Ok(res_vec)
}

/// get allele frequency from a list of vcf files
///
/// Only biallic SNPs are considerred; others are silently ignored
///
/// the returned results is a hashmap, with key being chromosome name,
/// and the value for each chromosme is a vector of 2-tuple (chromosomal bp position and the allele frequency)
pub fn get_afreq_from_vcf(vcf_files: &[String]) -> Result<HashMap<String, Vec<(u32, f32)>>> {
    let mut res = HashMap::<String, Vec<(u32, f32)>>::new();

    for vcf_file in vcf_files.iter() {
        let mut last_chr_id = u32::MAX;
        let mut v_place_holder = vec![];
        let mut last_v = &mut v_place_holder;
        let mut reader = bcf::Reader::from_path(vcf_file).context(HtslibSnafu {})?;
        let header = reader.header().clone();
        let mut rec = reader.empty_record();
        while reader.read(&mut rec).is_some() {
            // only consider biallelic SNPs
            if !((rec.allele_count() == 2) && unsafe { bcf_is_snp(rec.inner) != 0 }) {
                continue;
            }
            let mut ref_alt_cnt = [0u32; 2];
            for gt_sample in rec.format(b"GT").integer().context(HtslibSnafu {})?.iter() {
                for gta in *gt_sample {
                    let a: GenotypeAllele = (*gta).into();
                    if let Some(i) = a.index() {
                        ref_alt_cnt[i as usize] += 1
                    }
                }
            }

            let total = ref_alt_cnt[0] + ref_alt_cnt[1];
            if total == 0 {
                continue;
            }
            let afreq = ref_alt_cnt[1] as f32 / (total as f32);

            let pos = rec.pos() as u32;
            let chrid = rec.rid().context(VcfMissingRidSnafu {})?;
            if chrid == last_chr_id {
                last_v.push((pos, afreq));
            } else {
                let last_chr_name =
                    std::str::from_utf8(header.rid2name(chrid).context(HtslibSnafu {})?)
                        .context(Utf8Snafu {})?;
                println!("chrname: {last_chr_name}");
                res.entry(last_chr_name.to_owned())
                    .and_modify(|v| v.push((pos, afreq)))
                    .or_insert(vec![(pos, afreq)]);
                last_chr_id = chrid;
                // unwrap_or the default will never be used
                // this is only to make compiler happy
                last_v = res.get_mut(last_chr_name).unwrap_or(&mut v_place_holder);
            }
        }
    }
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::GenomeInfo;
    use ahash::HashMap;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Helper function to create a simple test VCF file
    fn create_test_vcf() -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");

        let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1	1/1	0/1
chr1	200	.	G	C	60	PASS	.	GT	0/0	0/0	0/1	1/1
chr1	300	.	C	G	60	PASS	.	GT	1/1	1/1	1/1	1/1
"#;

        file.write_all(vcf_content.as_bytes())
            .expect("Failed to write VCF content");
        file
    }

    // Helper function to create VCF with non-biallelic variants
    fn create_non_biallelic_vcf() -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");

        let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T,G	60	PASS	.	GT	0/1	1/2
chr1	200	.	G	C	60	PASS	.	GT	0/0	0/1
"#;

        file.write_all(vcf_content.as_bytes())
            .expect("Failed to write VCF content");
        file
    }

    // Helper function to create VCF with indels
    fn create_indel_vcf() -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");

        let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	ATG	A	60	PASS	.	GT	0/0	0/1
chr1	200	.	G	C	60	PASS	.	GT	0/0	0/1
"#;

        file.write_all(vcf_content.as_bytes())
            .expect("Failed to write VCF content");
        file
    }

    // Helper function to create VCF with missing genotypes
    fn create_missing_gt_vcf() -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");

        let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
chr1	100	.	A	T	60	PASS	.	GT	./.	0/1	1/1
chr1	200	.	G	C	60	PASS	.	GT	0/0	./.	./.
"#;

        file.write_all(vcf_content.as_bytes())
            .expect("Failed to write VCF content");
        file
    }

    // Helper function to create multi-chromosome VCF
    fn create_multi_chr_vcf() -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");

        let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##contig=<ID=chr2,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/1
chr1	200	.	G	C	60	PASS	.	GT	0/1	1/1
chr2	150	.	T	A	60	PASS	.	GT	0/0	0/1
chr2	250	.	C	G	60	PASS	.	GT	1/1	1/1
"#;

        file.write_all(vcf_content.as_bytes())
            .expect("Failed to write VCF content");
        file
    }

    #[test]
    fn test_get_afreq_from_vcf_basic() {
        let vcf_file = create_test_vcf();
        let vcf_files = vec![vcf_file.path().to_string_lossy().to_string()];

        let result = get_afreq_from_vcf(&vcf_files).expect("Should successfully parse VCF");

        assert_eq!(result.len(), 1); // One chromosome
        assert!(result.contains_key("chr1"));

        let chr1_data = &result["chr1"];
        assert_eq!(chr1_data.len(), 3); // Three SNPs

        // Position 100: 0/0, 0/1, 1/1, 0/1 = 4 alt alleles out of 8 total = 0.5
        assert_eq!(chr1_data[0].0, 99); // BCF is 0-based, VCF 1-based
        assert!((chr1_data[0].1 - 0.5).abs() < 0.001);

        // Position 200: 0/0, 0/0, 0/1, 1/1 = 3 alt alleles out of 8 total = 0.375
        assert_eq!(chr1_data[1].0, 199);
        assert!((chr1_data[1].1 - 0.375).abs() < 0.001);

        // Position 300: all 1/1 = 0 ref, 8 alt = 1.0
        assert_eq!(chr1_data[2].0, 299);
        assert!((chr1_data[2].1 - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_get_afreq_from_vcf_multiple_files() {
        let vcf_file1 = create_test_vcf();
        let vcf_file2 = create_multi_chr_vcf();

        let vcf_files = vec![
            vcf_file1.path().to_string_lossy().to_string(),
            vcf_file2.path().to_string_lossy().to_string(),
        ];

        let result =
            get_afreq_from_vcf(&vcf_files).expect("Should successfully parse multiple VCFs");

        assert!(result.contains_key("chr1"));
        assert!(result.contains_key("chr2"));

        // Should have combined data from both files for chr1
        let chr1_data = &result["chr1"];
        assert!(!chr1_data.is_empty());
    }

    #[test]
    fn test_get_afreq_from_vcf_filters_non_biallelic() {
        let vcf_file = create_non_biallelic_vcf();
        let vcf_files = vec![vcf_file.path().to_string_lossy().to_string()];

        let result = get_afreq_from_vcf(&vcf_files).expect("Should successfully parse VCF");

        assert_eq!(result.len(), 1);
        let chr1_data = &result["chr1"];
        assert_eq!(chr1_data.len(), 1); // Only the biallelic SNP should be included
        assert_eq!(chr1_data[0].0, 199); // Position 200 (0-based: 199)
    }

    #[test]
    fn test_get_afreq_from_vcf_filters_indels() {
        let vcf_file = create_indel_vcf();
        let vcf_files = vec![vcf_file.path().to_string_lossy().to_string()];

        let result = get_afreq_from_vcf(&vcf_files).expect("Should successfully parse VCF");

        assert_eq!(result.len(), 1);
        let chr1_data = &result["chr1"];
        assert_eq!(chr1_data.len(), 1); // Only the SNP should be included
        assert_eq!(chr1_data[0].0, 199); // Position 200 (0-based: 199)
    }

    #[test]
    fn test_get_afreq_from_vcf_handles_missing_genotypes() {
        let vcf_file = create_missing_gt_vcf();
        let vcf_files = vec![vcf_file.path().to_string_lossy().to_string()];

        let result = get_afreq_from_vcf(&vcf_files).expect("Should successfully parse VCF");

        assert_eq!(result.len(), 1);
        let chr1_data = &result["chr1"];

        // Both positions should be included because both have at least one valid genotype
        assert_eq!(chr1_data.len(), 2); // Two sites with valid genotypes

        // Position 100: ./., 0/1, 1/1 = 0 ref, 3 alt out of 4 total = 0.75
        assert_eq!(chr1_data[0].0, 99);
        assert!((chr1_data[0].1 - 0.75).abs() < 0.001);

        // Position 200: 0/0, ./., ./. = 2 ref, 0 alt out of 2 total = 0.0
        assert_eq!(chr1_data[1].0, 199);
        assert!((chr1_data[1].1 - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_get_afreq_from_vcf_empty_file_list() {
        let vcf_files: Vec<String> = vec![];
        let result = get_afreq_from_vcf(&vcf_files).expect("Should handle empty file list");
        assert!(result.is_empty());
    }

    #[test]
    fn test_get_afreq_from_vcf_nonexistent_file() {
        let vcf_files = vec!["nonexistent_file.vcf".to_string()];
        let result = get_afreq_from_vcf(&vcf_files);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), Error::HtslibError { .. }));
    }

    #[test]
    fn test_get_afreq_from_vcf_genome_wide_basic() {
        let vcf_file = create_multi_chr_vcf();
        let vcf_files = vec![vcf_file.path().to_string_lossy().to_string()];

        // Create a simple genome info for testing
        let mut idx = HashMap::new();
        idx.insert("chr1".to_string(), 0);
        idx.insert("chr2".to_string(), 1);
        let ginfo = GenomeInfo::new_from_parts(
            "test".to_string(),
            vec![1000, 1000],
            vec!["chr1".to_string(), "chr2".to_string()],
            idx,
            vec![0, 1000], // genome-wide starts: chr1 starts at 0, chr2 starts at 1000
        );

        let result = get_afreq_from_vcf_genome_wide(&vcf_files, &ginfo)
            .expect("Should successfully convert to genome-wide coordinates");

        assert!(!result.is_empty());

        // Results should be sorted by genome-wide position
        for i in 1..result.len() {
            assert!(
                result[i - 1].0 <= result[i].0,
                "Results should be sorted by genome-wide position"
            );
        }

        // Check that positions are correctly converted
        // chr1 positions should be < 1000 (first chromosome size)
        // chr2 positions should be >= 1000
        let chr1_positions: Vec<_> = result.iter().filter(|(pos, _)| *pos < 1000).collect();
        let chr2_positions: Vec<_> = result.iter().filter(|(pos, _)| *pos >= 1000).collect();

        assert!(!chr1_positions.is_empty());
        assert!(!chr2_positions.is_empty());
    }

    #[test]
    fn test_get_afreq_from_vcf_genome_wide_empty_input() {
        let vcf_files: Vec<String> = vec![];
        let mut idx = HashMap::new();
        idx.insert("chr1".to_string(), 0);
        let ginfo = GenomeInfo::new_from_parts(
            "test".to_string(),
            vec![1000],
            vec!["chr1".to_string()],
            idx,
            vec![0],
        );

        let result =
            get_afreq_from_vcf_genome_wide(&vcf_files, &ginfo).expect("Should handle empty input");
        assert!(result.is_empty());
    }

    #[test]
    fn test_get_afreq_from_vcf_genome_wide_propagates_errors() {
        let vcf_files = vec!["nonexistent_file.vcf".to_string()];
        let mut idx = HashMap::new();
        idx.insert("chr1".to_string(), 0);
        let ginfo = GenomeInfo::new_from_parts(
            "test".to_string(),
            vec![1000],
            vec!["chr1".to_string()],
            idx,
            vec![0],
        );

        let result = get_afreq_from_vcf_genome_wide(&vcf_files, &ginfo);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), Error::HtslibError { .. }));
    }

    #[test]
    fn test_error_types_can_be_created() {
        // Test that we can create each error type
        let invalid_bytes = vec![0xFF, 0xFE]; // Invalid UTF-8 sequence
        let utf8_error = std::str::from_utf8(&invalid_bytes).unwrap_err();
        let _utf8_snafu_error = Error::Utf8Error {
            source: utf8_error,
            backtrace: None,
        };

        let _vcf_missing_rid_error = Error::VcfMissingRid;

        // HtslibError would require creating an actual htslib error, which is complex
        // so we'll just verify it can be pattern matched
        let test_error = Error::VcfMissingRid;
        match test_error {
            Error::VcfMissingRid => (),
            Error::HtslibError { .. } => (),
            Error::Utf8Error { .. } => (),
        }
    }

    #[test]
    fn test_result_type_alias() {
        // Test that our Result type alias works correctly
        fn returns_result() -> Result<i32> {
            Ok(42)
        }

        let result = returns_result();
        assert_eq!(result.unwrap(), 42);
    }

    #[cfg(test)]
    mod edge_cases {
        use super::*;

        #[test]
        fn test_vcf_with_only_reference_calls() {
            let mut file = NamedTempFile::new().expect("Failed to create temp file");
            let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	0/0	0/0
"#;
            file.write_all(vcf_content.as_bytes())
                .expect("Failed to write VCF content");

            let vcf_files = vec![file.path().to_string_lossy().to_string()];
            let result =
                get_afreq_from_vcf(&vcf_files).expect("Should handle reference-only calls");

            let chr1_data = &result["chr1"];
            assert_eq!(chr1_data.len(), 1);
            assert!((chr1_data[0].1 - 0.0).abs() < 0.001); // Should be 0.0 allele frequency
        }

        #[test]
        fn test_vcf_with_only_alternate_calls() {
            let mut file = NamedTempFile::new().expect("Failed to create temp file");
            let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	1/1	1/1
"#;
            file.write_all(vcf_content.as_bytes())
                .expect("Failed to write VCF content");

            let vcf_files = vec![file.path().to_string_lossy().to_string()];
            let result =
                get_afreq_from_vcf(&vcf_files).expect("Should handle alternate-only calls");

            let chr1_data = &result["chr1"];
            assert_eq!(chr1_data.len(), 1);
            assert!((chr1_data[0].1 - 1.0).abs() < 0.001); // Should be 1.0 allele frequency
        }

        #[test]
        fn test_vcf_with_all_missing_genotypes_skipped() {
            let mut file = NamedTempFile::new().expect("Failed to create temp file");
            let vcf_content = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1500000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2
chr1	100	.	A	T	60	PASS	.	GT	./.	./.
chr1	200	.	G	C	60	PASS	.	GT	0/1	1/1
"#;
            file.write_all(vcf_content.as_bytes())
                .expect("Failed to write VCF content");

            let vcf_files = vec![file.path().to_string_lossy().to_string()];
            let result = get_afreq_from_vcf(&vcf_files).expect("Should skip all-missing sites");

            let chr1_data = &result["chr1"];
            assert_eq!(chr1_data.len(), 1); // Only the second site should be included
            assert_eq!(chr1_data[0].0, 199); // Position 200 (0-based: 199)
        }
    }

    #[cfg(test)]
    mod integration_with_testdata {
        use super::*;
        use std::path::PathBuf;

        fn get_testdata_path() -> PathBuf {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../testdata/dir001")
        }

        #[test]
        fn test_with_real_bcf_files() {
            let testdata_path = get_testdata_path();
            if !testdata_path.exists() {
                eprintln!("Skipping test: testdata not found at {testdata_path:?}");
                return;
            }

            let bcf_file = testdata_path.join("bcf/sel_chr1.bcf");
            if !bcf_file.exists() {
                eprintln!("Skipping test: BCF file not found");
                return;
            }

            let vcf_files = vec![bcf_file.to_string_lossy().to_string()];
            let result = get_afreq_from_vcf(&vcf_files);

            match result {
                Ok(data) => {
                    assert!(!data.is_empty(), "Should have chromosome data");
                    for (chrom, positions) in &data {
                        println!("Chromosome {}: {} positions", chrom, positions.len());
                        assert!(
                            !positions.is_empty(),
                            "Should have position data for {chrom}",
                        );

                        // Verify positions are sorted
                        for i in 1..positions.len() {
                            assert!(
                                positions[i - 1].0 <= positions[i].0,
                                "Positions should be sorted in chromosome {chrom}",
                            );
                        }

                        // Verify allele frequencies are in valid range
                        for (pos, freq) in positions {
                            assert!(
                                *freq >= 0.0 && *freq <= 1.0,
                                "Invalid allele frequency {freq} at position {pos}",
                            );
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Test failed with error: {e:?}");
                    panic!("Failed to process real BCF file");
                }
            }
        }

        #[test]
        fn test_with_genome_info_from_testdata() {
            let testdata_path = get_testdata_path();
            if !testdata_path.exists() {
                eprintln!("Skipping test: testdata not found");
                return;
            }

            let genome_toml = testdata_path.join("genome.toml");
            let bcf_file = testdata_path.join("bcf/sel_chr1.bcf");

            if !genome_toml.exists() || !bcf_file.exists() {
                eprintln!("Skipping test: required files not found");
                return;
            }

            let ginfo = match GenomeInfo::from_toml_file(genome_toml.to_str().unwrap()) {
                Ok(g) => g,
                Err(_) => {
                    eprintln!("Skipping test: could not load genome info");
                    return;
                }
            };

            let vcf_files = vec![bcf_file.to_string_lossy().to_string()];
            let result = get_afreq_from_vcf_genome_wide(&vcf_files, &ginfo);

            match result {
                Ok(data) => {
                    assert!(!data.is_empty(), "Should have genome-wide position data");

                    // Verify sorting by genome-wide position
                    for i in 1..data.len() {
                        assert!(
                            data[i - 1].0 <= data[i].0,
                            "Genome-wide positions should be sorted"
                        );
                    }

                    // Verify allele frequencies
                    for (gw_pos, freq) in &data {
                        assert!(
                            *freq >= 0.0 && *freq <= 1.0,
                            "Invalid allele frequency {freq} at genome-wide position {gw_pos}",
                        );
                    }

                    println!("Successfully processed {} variants", data.len());
                }
                Err(e) => {
                    eprintln!("Test failed with error: {e:?}");
                    panic!("Failed genome-wide conversion test");
                }
            }
        }
    }
}
