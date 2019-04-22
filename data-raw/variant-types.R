# VEP annotation as coded by vcf2maf (https://github.com/mskcc/vcf2maf)

coding_mutations = c('Frame_Shift_Del',
                     'Frame_Shift_Ins',
                     'In_Frame_Del',
                     'In_Frame_Ins',
                     'Missense_Mutation',
                     'Nonsense_Mutation',
                     'Nonstop_Mutation',
                     'Splice_Site',
                     'Targeted_Region',
                     'Translation_Start_Site')

truncating_mutations = c('Nonsense_Mutation',
                         'Frame_Shift_Ins',
                         'Frame_Shift_Del',
                         'Splice_Site',
                         'Nonstop_Mutation')

usethis::use_data(coding_mutations, truncating_mutations)
