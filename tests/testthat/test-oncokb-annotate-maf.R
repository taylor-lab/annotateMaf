context('oncokb_annotate_maf')
suppressPackageStartupMessages({
    library(tidyverse)
})
# test_that('warns if missing columns in input', {
#     
# })

test_that('returns valid output', {
    
    test_maf = tribble(~Hugo_Symbol, ~HGVSp_Short, ~HGVSc, ~Variant_Classification, ~Protein_position,
                       'KRAS', 'p.G12D', NULL, 'Missense_Mutation', '12/',
                       'TTN', 'p.S10275V', NULL, 'Missense_Mutation', '10275/109224')
    test_output = oncokb_annotate_maf(test_maf)
    
    expect_is(test_output, 'data.frame')
    expect_equal(ncol(test_output), ncol(test_maf) + 6)
    expect_identical(
        c(names(test_maf), 'cancer_type', 'oncogenic', 'oncokb_level',
                           'oncokb_resistance_level', 'oncokb_drugs', 'oncokb_version'),
                    names(test_output)
        )
    expect_match(test_output$oncogenic[which(test_output$Hugo_Symbol == 'KRAS' &
                                                 test_output$HGVSp_Short == 'p.G12D')], 'Oncogenic')
})
