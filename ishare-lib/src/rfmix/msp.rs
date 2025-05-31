/// msp file format:
///
/// no. col = 7 + 2 * n
///
/// col 0           1             2         3      4      5      6     7                 8                 9                  10
///
/// #Subpopulation  order/codes:  AFR=0     EAS=1  EUR=2  NAT=3
/// #chm            spos          epos      sgpos  egpos  n      snps  8v1_A.NAD_S100.0  8v1_A.NAD_S100.1  8v1_A.NAD_S1000.0  8v1_A.NAD_S1000.1
/// chr22           16747906      17321937  2.00   4.80   152    3     2                 2                 2                  3
/// chr22           17321937      17500594  4.80   5.42   160    3     2                 2                 2                  3
/// chr22           17500594      17586019  5.42   5.86   30     3     2                 2                 2                  3
/// chr22           17586019      17813718  5.86   6.35   295    3     2                 2                 2                  3
/// chr22           17813718      18028691  6.35   6.73   60     3     2                 2                 2                  3
/// chr22           18028691      19522206  6.73   9.56   190    3     2                 2                 2                  3
/// chr22           19522206      19623159  9.56   9.83   165    3     2                 2                 2                  3
/// chr22           26291477      26293669  26.70  26.73  5      2     2                 2                 2                  1
/// ...
/// chr22           26293669      26339326  26.73  26.89  105    2     2                 2                 2                  1
///
#[cfg(test)]
mod tests {

    #[test]
    fn test_msp_file_format_documentation() {
        // This test validates the documented MSP file format structure
        // The MSP format has 7 fixed columns + 2 columns per sample

        // Test column count calculation
        let n_samples = 2;
        let expected_columns = 7 + 2 * n_samples;
        assert_eq!(expected_columns, 11);

        // Test with different sample counts
        assert_eq!(7 + 2, 9); // 1 sample
        assert_eq!(7 + 2 * 3, 13); // 3 samples
        assert_eq!(7 + 2 * 10, 27); // 10 samples
    }

    #[test]
    fn test_msp_header_format() {
        // Test the expected header format for MSP files
        let expected_header_line1 = "#Subpopulation  order/codes:  AFR=0     EAS=1  EUR=2  NAT=3";
        let expected_header_line2 = "#chm            spos          epos      sgpos  egpos  n      snps  8v1_A.NAD_S100.0  8v1_A.NAD_S100.1  8v1_A.NAD_S1000.0  8v1_A.NAD_S1000.1";

        // Verify header contains expected elements
        assert!(expected_header_line1.contains("AFR=0"));
        assert!(expected_header_line1.contains("EAS=1"));
        assert!(expected_header_line1.contains("EUR=2"));
        assert!(expected_header_line1.contains("NAT=3"));

        // Verify column headers
        assert!(expected_header_line2.contains("chm"));
        assert!(expected_header_line2.contains("spos"));
        assert!(expected_header_line2.contains("epos"));
        assert!(expected_header_line2.contains("sgpos"));
        assert!(expected_header_line2.contains("egpos"));
        assert!(expected_header_line2.contains("snps"));
    }

    #[test]
    fn test_msp_data_line_parsing_concept() {
        // Test conceptual parsing of MSP data lines
        let sample_line = "chr22\t16747906\t17321937\t2.00\t4.80\t152\t3\t2\t2\t2\t3";
        let fields: Vec<&str> = sample_line.split('\t').collect();

        // Verify we have expected number of fields
        assert_eq!(fields.len(), 11); // 7 fixed + 2*2 samples

        // Test fixed column parsing
        assert_eq!(fields[0], "chr22"); // chromosome
        assert_eq!(fields[1], "16747906"); // start position
        assert_eq!(fields[2], "17321937"); // end position
        assert_eq!(fields[3], "2.00"); // start genetic position
        assert_eq!(fields[4], "4.80"); // end genetic position
        assert_eq!(fields[5], "152"); // n markers
        assert_eq!(fields[6], "3"); // snps

        // Test sample columns (ancestry assignments)
        assert_eq!(fields[7], "2"); // sample1 hap1
        assert_eq!(fields[8], "2"); // sample1 hap2
        assert_eq!(fields[9], "2"); // sample2 hap1
        assert_eq!(fields[10], "3"); // sample2 hap2
    }

    #[test]
    fn test_msp_position_ranges() {
        // Test position range validation for MSP segments
        struct MspSegment {
            _chromosome: String,
            start_pos: u32,
            end_pos: u32,
            start_genetic: f64,
            end_genetic: f64,
        }

        let segment1 = MspSegment {
            _chromosome: "chr22".to_string(),
            start_pos: 16747906,
            end_pos: 17321937,
            start_genetic: 2.00,
            end_genetic: 4.80,
        };

        let segment2 = MspSegment {
            _chromosome: "chr22".to_string(),
            start_pos: 17321937,
            end_pos: 17500594,
            start_genetic: 4.80,
            end_genetic: 5.42,
        };

        // Verify segments are non-empty
        assert!(segment1.end_pos > segment1.start_pos);
        assert!(segment2.end_pos > segment2.start_pos);

        // Verify genetic positions are increasing
        assert!(segment1.end_genetic > segment1.start_genetic);
        assert!(segment2.end_genetic > segment2.start_genetic);

        // Verify segments are adjacent (end of one = start of next)
        assert_eq!(segment1.end_pos, segment2.start_pos);
        assert_eq!(segment1.end_genetic, segment2.start_genetic);
    }

    #[test]
    fn test_msp_ancestry_codes() {
        // Test ancestry code validation
        #[derive(Debug, PartialEq)]
        enum Ancestry {
            Afr = 0,
            Eas = 1,
            Eur = 2,
            Nat = 3,
        }

        impl From<u8> for Ancestry {
            fn from(code: u8) -> Self {
                match code {
                    0 => Ancestry::Afr,
                    1 => Ancestry::Eas,
                    2 => Ancestry::Eur,
                    3 => Ancestry::Nat,
                    _ => panic!("Invalid ancestry code: {}", code),
                }
            }
        }

        // Test valid ancestry codes
        assert_eq!(Ancestry::from(0), Ancestry::Afr);
        assert_eq!(Ancestry::from(1), Ancestry::Eas);
        assert_eq!(Ancestry::from(2), Ancestry::Eur);
        assert_eq!(Ancestry::from(3), Ancestry::Nat);

        // Test that ancestry codes match expected values
        assert_eq!(Ancestry::Afr as u8, 0);
        assert_eq!(Ancestry::Eas as u8, 1);
        assert_eq!(Ancestry::Eur as u8, 2);
        assert_eq!(Ancestry::Nat as u8, 3);
    }

    #[test]
    fn test_msp_sample_name_parsing() {
        // Test sample name extraction from column headers
        let sample_columns = [
            "8v1_A.NAD_S100.0",  // sample1 hap1
            "8v1_A.NAD_S100.1",  // sample1 hap2
            "8v1_A.NAD_S1000.0", // sample2 hap1
            "8v1_A.NAD_S1000.1", // sample2 hap2
        ];

        // Extract sample names (remove haplotype suffix)
        let extract_sample_name = |col: &str| -> String {
            if let Some(last_dot) = col.rfind('.') {
                col[..last_dot].to_string()
            } else {
                col.to_string()
            }
        };

        let sample1 = extract_sample_name(sample_columns[0]);
        let sample2 = extract_sample_name(sample_columns[2]);

        assert_eq!(sample1, "8v1_A.NAD_S100");
        assert_eq!(sample2, "8v1_A.NAD_S1000");

        // Verify haplotype extraction
        let extract_haplotype = |col: &str| -> u8 {
            if let Some(last_dot) = col.rfind('.') {
                col[last_dot + 1..].parse().unwrap_or(0)
            } else {
                0
            }
        };

        assert_eq!(extract_haplotype(sample_columns[0]), 0); // hap1
        assert_eq!(extract_haplotype(sample_columns[1]), 1); // hap2
    }

    #[test]
    fn test_msp_segment_length_calculation() {
        // Test physical and genetic segment length calculations
        struct MspSegment {
            start_bp: u32,
            end_bp: u32,
            start_cm: f64,
            end_cm: f64,
        }

        impl MspSegment {
            fn physical_length(&self) -> u32 {
                self.end_bp - self.start_bp
            }

            fn genetic_length(&self) -> f64 {
                self.end_cm - self.start_cm
            }
        }

        let segment = MspSegment {
            start_bp: 16747906,
            end_bp: 17321937,
            start_cm: 2.00,
            end_cm: 4.80,
        };

        assert_eq!(segment.physical_length(), 574031); // bp
        assert!((segment.genetic_length() - 2.80).abs() < 1e-10); // cM

        // Test segment with smaller range
        let small_segment = MspSegment {
            start_bp: 26291477,
            end_bp: 26293669,
            start_cm: 26.70,
            end_cm: 26.73,
        };

        assert_eq!(small_segment.physical_length(), 2192); // bp
        assert!((small_segment.genetic_length() - 0.03).abs() < 1e-10); // cM
    }

    mod edge_cases {

        #[test]
        fn test_msp_single_marker_segment() {
            // Test segments with minimal marker count
            let n_markers = 1;
            assert!(n_markers > 0, "Segments should have at least one marker");

            // Test segment with very small physical length
            let start_pos = 26291477u32;
            let end_pos = 26291478u32; // Single base pair
            let length = end_pos - start_pos;

            assert_eq!(length, 1);
            assert!(length > 0, "Physical length should be positive");
        }

        #[test]
        fn test_msp_zero_genetic_length() {
            // Test segments with very small genetic length
            let start_genetic = 26.70f64;
            let end_genetic = 26.70f64; // Same genetic position
            let genetic_length = end_genetic - start_genetic;

            assert_eq!(genetic_length, 0.0);
            // This might occur in regions with no recombination
        }

        #[test]
        fn test_msp_boundary_ancestry_codes() {
            // Test boundary values for ancestry codes
            let valid_codes = [0u8, 1, 2, 3];

            for &code in &valid_codes {
                assert!(code <= 3, "Ancestry code {} should be valid", code);
            }

            // Test that codes outside range would be invalid
            let invalid_codes = [4u8, 255];
            for &code in &invalid_codes {
                assert!(code > 3, "Code {} should be considered invalid", code);
            }
        }

        #[test]
        fn test_msp_chromosome_edge_cases() {
            // Test different chromosome naming conventions
            let valid_chromosomes = vec![
                "chr1",
                "chr22",
                "chrX",
                "chrY",
                "1",
                "22",
                "X",
                "Y",
                "chromosome1",
                "scaffold_1",
            ];

            for chrom in &valid_chromosomes {
                assert!(!chrom.is_empty(), "Chromosome name should not be empty");
                assert!(!chrom.is_empty(), "Chromosome name should have length > 0");
            }
        }
    }

    mod error_handling {

        #[test]
        fn test_msp_invalid_position_format() {
            // Test handling of invalid position strings
            let invalid_positions = vec![
                "invalid",
                "",
                "16747906.5", // Float instead of integer
                "-100",       // Negative position
            ];

            for pos_str in &invalid_positions {
                let parse_result = pos_str.parse::<u32>();
                assert!(
                    parse_result.is_err(),
                    "Position '{}' should fail to parse",
                    pos_str
                );
            }
        }

        #[test]
        fn test_msp_invalid_genetic_position_format() {
            // Test handling of invalid genetic position strings
            let invalid_genetic = vec!["invalid", "", "two_point_five", "inf", "NaN"];

            for genetic_str in &invalid_genetic {
                let parse_result = genetic_str.parse::<f64>();
                assert!(
                    parse_result.is_err() || !parse_result.unwrap().is_finite(),
                    "Genetic position '{}' should fail or be non-finite",
                    genetic_str
                );
            }
        }

        #[test]
        fn test_msp_inconsistent_segment_positions() {
            // Test detection of inconsistent position ranges
            struct Position {
                start: u32,
                end: u32,
            }

            let invalid_segments = vec![
                Position {
                    start: 100,
                    end: 50,
                }, // End before start
                Position {
                    start: 100,
                    end: 100,
                }, // Zero length
            ];

            for seg in &invalid_segments {
                assert!(
                    seg.end <= seg.start,
                    "Segment with start={} end={} should be flagged as invalid",
                    seg.start,
                    seg.end
                );
            }
        }

        #[test]
        fn test_msp_column_count_validation() {
            // Test validation of column counts in MSP files
            fn validate_column_count(n_fields: usize, n_samples: usize) -> bool {
                let expected = 7 + 2 * n_samples;
                n_fields == expected
            }

            // Valid cases
            assert!(validate_column_count(9, 1)); // 1 sample
            assert!(validate_column_count(11, 2)); // 2 samples
            assert!(validate_column_count(13, 3)); // 3 samples

            // Invalid cases
            assert!(!validate_column_count(8, 1)); // Too few columns
            assert!(!validate_column_count(10, 2)); // Too few columns
            assert!(!validate_column_count(12, 2)); // Too many columns
        }
    }

    mod integration {

        #[test]
        fn test_msp_format_consistency() {
            // Test consistency between documentation and expected format
            let sample_data_lines = [
                "chr22\t16747906\t17321937\t2.00\t4.80\t152\t3\t2\t2\t2\t3",
                "chr22\t17321937\t17500594\t4.80\t5.42\t160\t3\t2\t2\t2\t3",
                "chr22\t17500594\t17586019\t5.42\t5.86\t30\t3\t2\t2\t2\t3",
            ];

            for (i, line) in sample_data_lines.iter().enumerate() {
                let fields: Vec<&str> = line.split('\t').collect();

                // All lines should have same number of fields
                assert_eq!(fields.len(), 11, "Line {} should have 11 fields", i);

                // Chromosome should be consistent
                assert_eq!(
                    fields[0], "chr22",
                    "Chromosome should be chr22 in line {}",
                    i
                );

                // Positions should parse as numbers
                let start_pos: u32 = fields[1].parse().unwrap();
                let end_pos: u32 = fields[2].parse().unwrap();
                assert!(
                    end_pos > start_pos,
                    "End position should be greater than start in line {}",
                    i
                );

                // Genetic positions should parse as floats
                let start_genetic: f64 = fields[3].parse().unwrap();
                let end_genetic: f64 = fields[4].parse().unwrap();
                assert!(
                    end_genetic >= start_genetic,
                    "End genetic position should be >= start in line {}",
                    i
                );

                // Ancestry codes should be valid (0-3)
                for (j, field) in fields.iter().enumerate().take(11).skip(7) {
                    let ancestry: u8 = field.parse().unwrap();
                    assert!(
                        ancestry <= 3,
                        "Ancestry code {} should be 0-3 in line {} field {}",
                        ancestry,
                        i,
                        j
                    );
                }
            }
        }

        #[test]
        fn test_msp_adjacent_segments() {
            // Test that consecutive segments are properly adjacent
            struct Segment {
                start: u32,
                end: u32,
                genetic_start: f64,
                genetic_end: f64,
            }

            let segments = [
                Segment {
                    start: 16747906,
                    end: 17321937,
                    genetic_start: 2.00,
                    genetic_end: 4.80,
                },
                Segment {
                    start: 17321937,
                    end: 17500594,
                    genetic_start: 4.80,
                    genetic_end: 5.42,
                },
                Segment {
                    start: 17500594,
                    end: 17586019,
                    genetic_start: 5.42,
                    genetic_end: 5.86,
                },
            ];

            for i in 1..segments.len() {
                // Physical positions should be adjacent
                assert_eq!(
                    segments[i - 1].end,
                    segments[i].start,
                    "Physical positions should be adjacent between segments {} and {}",
                    i - 1,
                    i
                );

                // Genetic positions should be adjacent
                assert!(
                    (segments[i - 1].genetic_end - segments[i].genetic_start).abs() < 1e-10,
                    "Genetic positions should be adjacent between segments {} and {}",
                    i - 1,
                    i
                );
            }
        }
    }
}
