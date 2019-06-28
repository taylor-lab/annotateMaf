#' Annotate somatic hotspot mutations
#'
#' Adds hotspot annotation to VEP-annotated MAF. Sources of default hotspots below.
#'
#' @param maf Input MAF.
#' @param hotspot_tbl Custom table of hotspots.
#'
#' @return Annotated MAF with columns \code{snv_hotspot}, \code{threeD_hotspot}, \code{indel_hotspot_type} and \code{Hotspot} indicating types of hotspots. Note that the \code{Hotspot} column does not includes 3D hotspots.
#'
#' @importFrom dplyr mutate rowwise case_when pull filter
#' @importFrom purrr discard map
#' @importFrom tidyr replace_na
#' @importFrom jsonlite fromJSON
#' @importFrom readr read_lines
#' @importFrom data.table %like%
#' @importFrom stringr str_c str_split str_extract str_replace
#'
#' @name hotspot_annotate_maf
NULL

load_gene_annotation = function() {
    read_lines('http://oncokb.org/api/v1/genes') %>%
        fromJSON()
}

#' @export
#' @rdname hotspot_annotate_maf
hotspot_annotate_maf = function(maf, hotspot_tbl = NULL) {
    
    if (!inherits(maf, 'data.frame')) stop('Input MAF must be a data frame, preferrable VEP annotated')
    if (!is.null(hotspot_tbl)) {
        if (!file.exists(hotspot_tbl)) {
            stop('Hotspots file does not exist', call. = F)
        } else {
            hotspots = hotspot_tbl
        }
    } else {
        hotspots = annotateMaf::hotspots
    }
    if (any(grepl('hotspot', tolower(names(maf))))) message('Hotspot columns in MAF might be overwritten, check names')
    
    # Read default hotspot lists if user does not supply
    gene_annotation = load_gene_annotation()

    # Function that deals with indel hotspots
    tag_indel_hotspot = function(gene, hgvsp_short, start, end, indel_length) {
        is_hotspot = 'none'
        gene_hotspots = dplyr::filter(hotspots, Gene == gene, indel_hotspot == T) %>% # checks if previously identified
            dplyr::pull(previous_mutations) %>%
            stringr::str_split(',') %>%
            unlist() %>%
            purrr::discard(. == '')
        start_res = as.numeric(stringr::str_extract(gene_hotspots, '[0-9]+(?=_)'))
        end_res = as.numeric(stringr::str_extract(gene_hotspots, '(?<=_[A-Z]{1})[0-9]+'))
        longest_variant = max(end_res - start_res, na.rm = TRUE)
        indel_hs = dplyr::filter(hotspots, Gene == gene, indel_hotspot == TRUE) %>% # checks if overlap with hotspot intervals
            dplyr::filter(Start - 1 <= (start + 2) & End >= (end - 2) & indel_length <= (longest_variant + 1)) # allow for wiggle room in both directions but not too long indel
        if (nrow(indel_hs) > 0 | stringr::str_replace(hgvsp_short, 'p.', '') %in% gene_hotspots) {
            is_hotspot = 'prior'
        } else if (end - start <= 5) { # if short hotspot overlapping known hotspot
            snv_hs = dplyr::filter(hotspots, Gene == gene & (indel_hotspot == FALSE | is.na(End))) %>% # includes indel hotspots that only cover one codon
                dplyr::filter(between(Pos, start, end))
            if (nrow(snv_hs) > 0) is_hotspot = 'novel'
        }
        return(is_hotspot)
    }

    tag_onp_hotspot = function(gene, type, trunc, start, end, mut_length) {
        is_hotspot = FALSE
        gene_hotspots = dplyr::filter(hotspots, Gene == gene, snv_hotspot == TRUE) %>% # checks if previously identified
            dplyr::select(tag, previous_mutations) %>%
            dplyr::mutate(previous_mutations = purrr::map(previous_mutations, ~unlist(stringr::str_split(., ','))))
        if (type == 'SNP' & trunc == FALSE) {
            is_hotspot = stringr::str_c(gene, start) %in% gene_hotspots$tag
        } else if (type == 'SNP' & trunc == TRUE) {
            is_hotspot = stringr::str_c(gene, start) %in% gene_hotspots$tag &
                gene %in% gene_annotation$hugoSymbol[gene_annotation$tsg == TRUE]
        } else if (mut_length > 3) {
            is_hotspot = FALSE
        } else if (any(c(stringr::str_c(gene, start), stringr::str_c(gene, end)) %in% gene_hotspots$tag) &
                   trunc == FALSE) {
            is_hotspot = TRUE
        } else if (any(c(stringr::str_c(gene, start), stringr::str_c(gene, end)) %in% gene_hotspots$tag) &
                   trunc == TRUE) {
            is_hotspot = gene %in% gene_annotation$hugoSymbol[gene_annotation$tsg == TRUE]
        }
        return(is_hotspot)
    }

    maf = dplyr::mutate(maf,
                        residue = stringr::str_extract(Protein_position, '^[0-9]+(?=/|-)'),
                        start_residue = residue,
                        end_residue = stringr::str_extract(Protein_position, '(?<=-)[0-9]+(?=/)'),
                        end_residue = ifelse(is.na(end_residue) & Variant_Classification %like% 'In_Frame',
                                             start_residue, end_residue)) %>%
        tidyr::replace_na(list(start_residue = 0, end_residue = 0)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            start_residue = as.numeric(start_residue),
            end_residue = as.numeric(end_residue),
            variant_length = dplyr::case_when(
                Variant_Type %in% c('SNP', 'DNP', 'TNP', 'ONP') ~ nchar(Reference_Allele)/3,
                Variant_Classification == 'In_Frame_Del' ~ nchar(Reference_Allele)/3,
                Variant_Classification == 'In_Frame_Ins' ~ nchar(Tumor_Seq_Allele2)/3
            ),
            snv_hotspot = ifelse(Variant_Type %like% 'NP$' &
                                     Variant_Classification %in% coding_mutations &
                                     !is.na(residue),
                                 tag_onp_hotspot(
                                     Hugo_Symbol,
                                     Variant_Type,
                                     Variant_Classification %in% truncating_mutations,
                                     as.numeric(start_residue),
                                     as.numeric(end_residue),
                                     variant_length
                                 ),
                                 FALSE),
            threeD_hotspot = ifelse(Variant_Type %in% c('SNP', 'DNP') &
                                        Variant_Classification %in% coding_mutations & !is.na(residue),
                                    stringr::str_c(Hugo_Symbol, residue) %in% hotspots$tag[hotspots$threeD_hotspot == T],
                                    FALSE),
            indel_hotspot_type = ifelse(Variant_Classification %like% 'In_Frame',
                                        tag_indel_hotspot(Hugo_Symbol,
                                                          HGVSp_Short,
                                                          as.numeric(start_residue),
                                                          as.numeric(end_residue),
                                                          variant_length),
                                        'none'),
            indel_hotspot = indel_hotspot_type != 'none',
            Hotspot = snv_hotspot == TRUE | indel_hotspot == TRUE)

    maf
}
