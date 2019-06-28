context('hotspot_annotate_maf')
suppressPackageStartupMessages({
    library(dplyr)
    library(tribble)
})

test_that('returns valid output', {
    
    test_maf = tribble(~Hugo_Symbol, ~Reference_Allele, ~Tumor_Seq_Allele2, ~HGVSp_Short, ~HGVSc,
                       ~Variant_Classification, ~Protein_position, ~Variant_Type,
                       'KRAS', 'G', 'A', 'p.G12D', NULL, 'Missense_Mutation', '12/', 'SNP',
                       'TTN', 'C', 'G', 'p.S10275V', NULL, 'Missense_Mutation', '10275/109224', 'SNP')
    test_output = hotspot_annotate_maf(test_maf)
    
    expect_is(test_output, 'data.frame')
    expect_equal(ncol(test_output), ncol(test_maf) + 9)
    expect_identical(
        c(names(test_maf), 'residue', 'start_residue', 'end_residue',
          'variant_length', 'snv_hotspot', 'threeD_hotspot', 'indel_hotspot_type', 'indel_hotspot', 'Hotspot'),
        names(test_output)
    )
    expect_true(test_output$Hotspot[which(test_output$Hugo_Symbol == 'KRAS' & test_output$HGVSp_Short == 'p.G12D')])
})
