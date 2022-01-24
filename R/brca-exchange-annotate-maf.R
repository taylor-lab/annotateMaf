#' Annotate BRCA1/2 variants
#'
#' Adds functional BRCA1 and BRCA2 variant annotation, using BRCA Exchange API (see URL). Annotation comes from ENIGMA and ClinVar.
#' Not that variant alleles in a MAF are in the column \code{Tumor_Seq_Allele2} by default.
#'
#' @param maf Input MAF.
#' 
#' @param gene Query gene.
#' @param start Variant start position, in hg19.
#' @param end Variant end position, in hg19.
#' @param ref Reference allele.
#' @param alt Alternate allele.
#' @param use_api If \code{TRUE} uses BRCA exchange API, otherwise static table.
#' 
#' @return Annotated MAF with columns \code{brca_exchange_enigma}, \code{brca_exchange_clinvar} variant annotation from ENIGMA and ClinVar, respectively, and \code{brca_exchange_id} indicating database variant ID.
#' 
#' @source \url{https://brcaexchange.org}
#'
#' @importFrom dplyr filter mutate select tibble case_when
#' @importFrom purrr map_dfr map_chr possibly transpose pmap set_names
#' @importFrom reticulate source_python
#' @importFrom stringr str_sub str_trim
#' @importFrom tibble tibble
#' 
#' @examples
#' \donttest{query_brca_exchange('BRCA1', 41276045, 41276046, 'CT', '-')}
#'   
#' @name brca_exchange_annotate_maf
NULL

#' @export
#' @rdname brca_exchange_annotate_maf
query_brca_exchange = function(gene, start, end, ref, alt, use_api = FALSE) {
    
    if (end < start) stop('End position cannot be smaller than start position.', call. = F)
    
    if (gene %in% c('BRCA1', 'BRCA2')) {
        
        if (use_api == TRUE) {
            qq = brca_query(gene, start - 1, end + 1)
        } else {
            chrom = ifelse(gene == 'BRCA1', 17, 13)
            qq = dplyr::filter(brca_exchange_variants, Chr == chrom, pyhgvs_Hg37_Start > start - 2, pyhgvs_Hg37_End < end + 2) %>%
                purrr::transpose()
        }
        
        if (length(qq) > 0) {
            purrr::map_dfr(qq, ~purrr::set_names(., c('gene',
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
                    variant_type = dplyr::case_when(
                        stringr::str_sub(reference_allele, 1, 1) == stringr::str_sub(alternate_allele, 1, 1) &
                            nchar(reference_allele) < nchar(alternate_allele) ~ 'insertion',
                        stringr::str_sub(reference_allele, 1, 1) == stringr::str_sub(alternate_allele, 1, 1) &
                            nchar(reference_allele) > nchar(alternate_allele) ~ 'deletion',
                        TRUE ~ 'snv'),
                    reference_allele = dplyr::case_when(
                        variant_type == 'insertion' ~ '-',
                        variant_type == 'deletion' ~ stringr::str_replace(reference_allele, '(^[A-Z]{1})', ''),
                        TRUE ~ reference_allele
                    ),
                    alternate_allele = dplyr::case_when(
                        variant_type == 'insertion' ~ stringr::str_replace(alternate_allele, '(^[A-Z]{1})', ''), 
                        variant_type == 'deletion' ~ '-',
                        TRUE ~ alternate_allele
                    ),
                    end_position = dplyr::case_when(
                        variant_type == 'snv' ~ start_position,
                        variant_type == 'insertion' ~ end_position + 1,
                        variant_type == 'deletion' ~ end_position),
                    start_position = dplyr::case_when(
                        variant_type == 'deletion' ~ start_position + 1,
                        TRUE ~ start_position)
                ) %>% 
                dplyr::filter(start_position == start,
                              end_position == end,
                              reference_allele == ref,
                              alternate_allele == alt) %>% 
                dplyr::mutate(brca_exchange_enigma = stringr::str_trim(str_extract(
                    pathogenicity, '[A-Za-z\\,\\ ]+(?=\\(ENIGMA\\))'), 'both'),
                    brca_exchange_clinvar = stringr::str_trim(str_extract(
                        pathogenicity, '[A-Za-z\\,\\ \\_]+(?=\\(ClinVar\\))'), 'both')) %>%
                dplyr::select(brca_exchange_id = id, brca_exchange_enigma, brca_exchange_clinvar)
        } else {
            tibble::tibble(brca_exchange_id = NA)
        }
    }
}

#' @export
#' @rdname brca_exchange_annotate_maf
brca_annotate_maf = function(maf, use_api = FALSE) {
    
    map_chr_possibly = purrr::possibly(purrr::map_chr, NA_character_)
    
    dplyr::mutate(maf, annot = purrr::pmap(list(Hugo_Symbol, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, use_api),
                                           query_brca_exchange)) %>%
        dplyr::mutate(brca_exchange_id = map_chr_possibly(annot, 'brca_exchange_id', .default = NA),
                      brca_exchange_enigma = map_chr_possibly(annot, 'brca_exchange_enigma', .default = NA),
                      brca_exchange_clinvar = map_chr_possibly(annot, 'brca_exchange_clinvar', .default = NA)) %>%
        dplyr::select(-annot)
}
