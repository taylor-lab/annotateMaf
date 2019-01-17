#' Hotspot
#'
#' Code used to generate lists of somatic hotspots.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{Gene}{Gene symbol.}
#'   \item{Residue}{Mutated protein residue}
#'   \item{Pos}{Transcript position.}
#'   \item{Start}{If indel, start position.}
#'   \item{End}{If indel, end position.}
#'   \item{indel_hotspot}{Boolean indicating if indel hotspot.}
#'   \item{snv_hotspot}{Boolean indicating if SNV hotspot.}
#'   \item{threeD_hotspot}{Boolean indicating if 3D hotspot.}
#'   \item{previous_mutations}{Comma-separate list of previously observed mutations at site.}
#'   \item{tag}{Tag used to match variant.}
#' }
#' 
#' import dplyr
#' import purrr
#' importFrom usethis use_data
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/26619011}
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/29247016}
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/28115009}
'hotspots'

# Generate hotspot list -------------------------------------------------------------------------------------------
# hs_from_local = c('/luna/work/hotspots/data/nbt-cd24k-FP-filtered-hotspots.txt',
#                   '/luna/work/hotspots/3d/hotspots.txt')
# 
# hs_from_cluster = c('/ifs/work/taylorlab/jonssonp/hotspots/data/nbt-cd24k-FP-filtered-hotspots.txt',
#                     '/ifs/work/taylorlab/jonssonp/hotspots/3d/hotspots.txt')
# 
# if (Sys.info()[['nodename']] == 'lski2423') {
#     hotspot_files =  hs_from_local
# } else {
#     hotspot_files = hs_from_cluster
# }
# 
# hotspots = map_dfr(hotspot_files, function(x) {
#     fread(x) %>%
#         mutate(source = x,
#                Type = if('Type' %in% names(.)) { Type } else { 'single residue' },
#                Class = if('Class' %in% names(.)) { Class } else { 'single residue' },
#                false_positive = if('false_positive' %in% names(.)) { false_positive } else { FALSE },
#                previous_mutations = if('previous_mutations' %in% names(.)) { previous_mutations } else { str_replace_all(Variants, ':[0-9]+\\|?', ',') })
# }) %>%
#     mutate(indel_hotspot = Type == 'in-frame indel',
#            indel_hotspot = ifelse(is.na(indel_hotspot), FALSE, indel_hotspot),
#            threeD_hotspot = ifelse(Class %in% c('Hotspot-linked', 'Cluster-exclusive'), TRUE, FALSE),
#            snv_hotspot = source %like% 'NBT|24k' & indel_hotspot == FALSE,
#            Pos = ifelse(indel_hotspot == F, str_extract(Residue, '(?<=[A-Z])[0-9]+'), NA),
#            Start = ifelse(indel_hotspot == T, str_extract(Residue, '^[0-9]+(?=-)'), NA),
#            End = ifelse(indel_hotspot == T, str_extract(Residue, '(?<=-)[0-9]+$'), NA),
#            Pos = ifelse(indel_hotspot == T & is.na(Start), Residue, Pos)) %>%
#     replace_na(list(false_positive = FALSE)) %>%
#     group_by(Gene, Residue, Pos, Start, End) %>%
#     dplyr::summarize(indel_hotspot = any(indel_hotspot == T, na.rm = T),
#                      snv_hotspot = any(snv_hotspot == T, na.rm = T),
#                      threeD_hotspot = any(threeD_hotspot == T, na.rm = T),
#                      previous_mutations = paste(c(unique(unlist(strsplit(previous_mutations, ',')))), collapse = ','),
#                      false_positive = all(false_positive == T, na.rm = T)) %>%
#     ungroup() %>%
#     mutate(tag = str_c(Gene, Pos),
#            Pos = as.numeric(Pos),
#            Start = as.numeric(Start),
#            End = as.numeric(End)) %>%
#     filter(false_positive == FALSE) %>% 
#     select(-false_positive)

# usethis::use_data(hotspots)
