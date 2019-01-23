#' Annotate BRCA1/2 variants
#'
#' Adds functional BRCA1 and BRCA2 variant annotation, using BRCA Exchange API (see URL). Annotation comes from ENIGMA and ClinVar.
#' Not that variant alleles in a MAF are in the column \code{Tumor_Seq_Allele2} by default.
#'
#' @param maf Input MAF.
#' 
#' @param gene Query gene.
#' @param start Variant start position.
#' @param end Variant end position.
#' @param ref Reference allele.
#' @param alt Alternate allele.
#' 
#' @return Annotated MAF with columns \code{brca_exchange_enigma}, \code{brca_exchange_clinvar} variant annotation from ENIGMA and ClinVar, respectively, and \code{brca_exchange_id} indicating database variant ID.
#' 
#' @source \url{https://brcaexchange.org}
#'
#' @import dplyr
#' @importFrom purrr map_dfr map_chr possibly
#' @importFrom reticulate source_python
#' 
#' @examples
#' \donttest{query_brca_exchange('BRCA1', 41276045, 41276046, 'CT', '-')}
#'   
#' @name brca_exchange_annotate_maf
NULL

#' @export
#' @rdname brca_exchange_annotate_maf
query_brca_exchange = function(gene, start, end, ref, alt) {
    
    if (end < start) stop('End position cannot be smaller than start position.', call. = F)
    
    if (gene %in% c('BRCA1', 'BRCA2')) {
        brca_query(gene, start-1, end+1) %>%
            map_dfr(., ~set_names(., c('gene',
                                       'chrom',
                                       'start_position',
                                       'end_position',
                                       'reference_allele',
                                       'alternate_allele',
                                       'id',
                                       'pathogenicity'))) %>% 
            dplyr::mutate(
                start_position = as.numeric(start_position),
                end_position = as.numeric(end_position),
                variant_type = case_when(
                    str_sub(reference_allele, 1, 1) == str_sub(alternate_allele, 1, 1) &
                        nchar(reference_allele) < nchar(alternate_allele) ~ 'insertion',
                    str_sub(reference_allele, 1, 1) == str_sub(alternate_allele, 1, 1) &
                        nchar(reference_allele) > nchar(alternate_allele) ~ 'deletion',
                    TRUE ~ 'snv'),
                reference_allele = case_when(
                    variant_type == 'insertion' ~ '-',
                    variant_type == 'deletion' ~ str_replace(reference_allele, '(^[A-Z]{1})', ''),
                    TRUE ~ reference_allele
                ),
                alternate_allele = case_when(
                    variant_type == 'insertion' ~ str_replace(alternate_allele, '(^[A-Z]{1})', ''), 
                    variant_type == 'deletion' ~ '-',
                    TRUE ~ alternate_allele
                ),
                end_position = case_when(
                    variant_type == 'snv' ~ start_position,
                    variant_type == 'insertion' ~ end_position + 1,
                    variant_type == 'deletion' ~ end_position),
                start_position = case_when(
                    variant_type == 'deletion' ~ start_position + 1,
                    TRUE ~ start_position)
            ) %>% 
            dplyr::filter(start_position == start, end_position == end, reference_allele == ref, alternate_allele == alt) %>% 
            mutate(brca_exchange_enigma = str_trim(str_extract(pathogenicity, '[A-Za-z\\,\\ ]+(?=\\(ENIGMA\\))'), 'both'),
                   brca_exchange_clinvar = str_trim(str_extract(pathogenicity, '[A-Za-z\\,\\ \\_]+(?=\\(ClinVar\\))'), 'both')) %>%
            select(brca_exchange_id = id, brca_exchange_enigma, brca_exchange_clinvar)
    }
}


#' @export
#' @rdname brca_exchange_annotate_maf
brca_annotate_maf = function(maf) {

    map_chr_possibly = possibly(map_chr, NA_character_)

    mutate(maf, annot = pmap(list(Hugo_Symbol, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2),
                             query_brca_exchange)) %>%
        mutate(brca_exchange_id = map_chr_possibly(annot, 'brca_exchange_id'),
               brca_exchange_enigma = map_chr_possibly(annot, 'brca_exchange_enigma'),
               brca_exchange_clinvar = map_chr_possibly(annot, 'brca_exchange_clinvar')) %>%
        select(-annot)
}
