#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#'
#' @export
hypGEM_to_excel <- function(
    hypeR_GEM_obj,
    file_path,
    cols = NULL,
    fdr_cutoff = 1.0,
    versioning = TRUE) {
  ## input checks
  stopifnot(.check_hypeR_GEM_obj(hypeR_GEM_obj))
  stopifnot(fdr_cutoff > 0 && fdr_cutoff <= 1.0)

  ## Generate excel file
  wb <- openxlsx::createWorkbook()

  ## A new sheet for each dataframe
  for ( nm in names(hypeR_GEM_obj) ) {
    sheet <- openxlsx::addWorksheet(wb, sheetName = nm)

    ## Extract hypGEM dataframe
    df <- hypeR_GEM_obj[[nm]]$data

    if (is.null(cols)) {
      cols <- seq_len(ncol(df))
    }
    df_flt <- df |>
      dplyr::select({{ cols }}) |>
      dplyr::filter(fdr <= fdr_cutoff)

    openxlsx::writeData(wb, sheet = nm, x = df_flt, colNames = TRUE, rowNames = FALSE)
  }
  # if ( versioning ) {
  #  x <- mapply(function(x) x$info, multihyp_obj$data)
  #  sheet <- openxlsx::addWorksheet(wb, sheetName = "versioning")
  #  openxlsx::writeData(wb, sheet = "versioning", x = x, colNames = ncol(x) > 1, rowNames = TRUE)
  # }
  suppressMessages(openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE))
}
