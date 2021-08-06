#' Get directory names produced by kallisto quant call.
#' 
#' @description
#' The directories created by kallisto quant (see my notebook entry
#'   on 20210623 for explanation of the exact call) contain the results
#'   we need in order to get the expression level percentiles.
#'
#' @param x The directory in which to begin the recursive search for
#'   kallisto output directories
#' @return Vector of directory names beginning with "SRR"
get_srr_dirs = function(x) {
  dirs = list.dirs(x)
  final_dirs = NA
  count = 1
  for (i in 1:length(dirs)) {
    dir = dirs[i]
    if (grepl("SRR", dir)) {
      final_dirs[count] = dir
      count = count + 1
    }
  }
  return(final_dirs)
}

#' Extract just the run id from a path
#' 
#' Get the SRR ID
#' @param x A path to split on '/' and grab just the SRA run accession.
#' @return The SRA accession for this path
grab_srr = function(x) {
  split_dir = str_split(x, '/', simplify=TRUE)
  col_num = dim(split_dir)[2]
  srr_name = split_dir[col_num]
  return(srr_name)
}

#' Convert vector of TPM values to percentiles
#' 
#' @description 
#' Written to be used on tibbles passed to this function by purrr::map,
#'   it grabs the vector of tpm values from the tibble, then
#'   calculates the percentile for each value.
#' 
#' @param x A tibble.
#' @return A vector of percentiles
calc_tpm_percentile = function(x) {
  tpm = x$tpm
  percs = NA
  for (i in 1:length(tpm)) {
    percs[i] = length(tpm[tpm <= tpm[i]]) / length(tpm) * 100
  }
  return(percs)
}

#' Get protein ID from the string present in NCBI's CDS fasta file
#' 
#' @description 
#' Splits up the fasta record name provided in NCBI's CDS fasta file
#'   and produces (and returns) a protein ID. This function was
#'   written to operate on tibbles passed by purrr::map.
#'   
#' @param x A tibble.
#' @return A vector of protein IDs
extract_prot_id = function(x) {
  ids = x$target_id
  splits = str_split(ids, '_', simplify=TRUE)[,c(4,5)]
  prot_ids = NA
  for (i in 1:nrow(splits)) {
    prot_ids[i] = paste(splits[i,], sep="_", collapse="_")
  }
  return(prot_ids)
}

#' Filter input tibble for GO terms of interest
#' 
#' @param in_tib a tibble.
#' @param uniprot_tib A tibble containing uniprot information.
#' @param YES_term GO term to retain in returned tibble.
#' @param NOT_terms A vector of GO terms to exclude from results.
#' @param perc_threshold Minimum percentile tpm to retain in returned tibble.
#' @return A filtered tibble
filter_NAPs = function(in_tib, uniprot_tib,
                       YES_term, NOT_terms,
                       perc_threshold=80) {
  
  tpms = in_tib %>%
    # join the crucial uniprot information to this tibble
    left_join(
      uniprot_tib %>% dplyr::select(
        protein,
        gene,
        locus_tag,
        GO_terms
      ),
      by="locus_tag"
    )
  
  binders_of_interest = tpms %>%
    dplyr::select(
      c(locus_tag, protein, gene, GO_terms, percentiles, condition, bug)
    ) %>%
    filter(
      grepl(YES_term, GO_terms),
      !grepl(paste(NOT_terms, collapse="|"), GO_terms)
    )
  
  potential_NAPs = binders_of_interest %>%
    group_by(locus_tag) %>%
    filter(
      any(percentiles >= perc_threshold)
    ) %>%
    ungroup() %>%
    arrange(desc(percentiles))
  return(potential_NAPs)
}

purr_widen = function(x, names, values) {
  names = enquo(names)
  values = enquo(values)
  new_df = x %>%
    tidyr::pivot_wider(names_from = !!names, values_from = !!values)
  return(new_df)
}