
# Pull list of BRCA variants from BRCA Exchange -------------------------------------------------------------------
# Thu Feb 14 09:44:53 2019 ------------------------------
# See https://github.com/BRCAChallenge/brca-exchange/issues/981

system('wget https://brcaexchange.org/backend/downloads/releases/current_release.tar.gz && tar -xzf current_release.tar.gz')
brca_exchange_variants = fread('output/release/built_with_change_types.tsv') %>%
    mutate(id = 'N/A') %>% 
    select(Gene_Symbol, Chr, pyhgvs_Hg37_Start, pyhgvs_Hg37_End, Ref, Alt, id, Pathogenicity_all)

system('rm current_release.tar.gz && rm -rf output')

usethis::use_data(brca_exchange_variants)
