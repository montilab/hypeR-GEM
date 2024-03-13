#' Writes the enrichment results out to an Excel workbook
#'

#' @param hypGEM_obj A hypeR-GEM object from enrichment analysis
#' @param file_path the path to save the output file
#' @param cols columns
#' @param fdr_cutoff fdr threshold
#' @param do_sort logical parameter
#' @param versioning logical parameter
#'
#' @return a excel file

#'
#' @import methods utils
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#'
#' @export
hypGEM_to_excel <- function(
    hypGEM_obj,
    file_path,
    cols = NULL,
    fdr_cutoff = 1.0,
    do_sort = FALSE,
    versioning = TRUE) {
  ## input checks
  stopifnot(.check_hypeR_GEM_obj(hypGEM_obj))
  stopifnot(fdr_cutoff > 0 && fdr_cutoff <= 1.0)

  ## Generate excel file
  wb <- openxlsx::createWorkbook()

  ## filter result tables
  if ( !is.null(cols) || fdr_cutoff < 1.0 )
    hypGEM_obj <- hypGEM_filter(hypGEM_obj, fdr_cutoff = fdr_cutoff, cols = cols, do_sort = do_sort)

  ## A new sheet for each dataframe
  for ( nm in names(hypGEM_obj) ) {
    sheet <- openxlsx::addWorksheet(wb, sheetName = nm)

    ## Extract hypGEM dataframe
    df <- hypGEM_obj[[nm]]$data

    openxlsx::writeData(wb, sheet = nm, x = df, colNames = TRUE, rowNames = FALSE)
  }
  # if ( versioning ) {
  #  x <- mapply(function(x) x$info, multihyp_obj$data)
  #  sheet <- openxlsx::addWorksheet(wb, sheetName = "versioning")
  #  openxlsx::writeData(wb, sheet = "versioning", x = x, colNames = ncol(x) > 1, rowNames = TRUE)
  # }
  suppressMessages(openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE))
}

#' Filter the enrichment results by fdr value, and optionally sub-selects result columns
#'
#' @param hypGEM_obj A hypeR-GEM object from enrichment analysis
#' @param cols columns
#' @param fdr_cutoff fdr threshold
#' @param do_sort logical parameter
#'
#' @return filtered data frame
#'
#' @importFrom rlang .data

#'
#' @export
hypGEM_filter <- function(
    hypGEM_obj,
    cols = NULL,
    fdr_cutoff = 1.0,
    do_sort = FALSE)
{
  ## input checks
  stopifnot(.check_hypeR_GEM_obj(hypGEM_obj))
  stopifnot(fdr_cutoff > 0 && fdr_cutoff <= 1.0)
  stopifnot( is.null(cols) || all(cols %in% hypGEM_obj[[1]]$data))

  ## nothing to do here
  if ( is.null(cols) && fdr_cutoff == 1 && !do_sort )
    return(hypGEM_obj)

  ## if no columns specified, select all
  if (is.null(cols)) {
    cols <- seq_len(ncol(hypGEM_obj[[1]]$data))
  }
  ## filter by fdr and by columns, and sort by fdr if required
  obj_flt <- hypGEM_obj |>
    purrr::map( \(ls) {
      ls$data <- ls$data |>
        dplyr::select({{ cols }}) |>
        dplyr::filter(.data$fdr <= fdr_cutoff)
      if (do_sort)
        ls$data <- ls$data |> dplyr::arrange(.data$pval)
      return(ls)
    })
  return(obj_flt)
}
