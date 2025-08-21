#' Format a string using placeholders
#'
#' @param string A an unformatted string with placeholders
#' @param ... Variables to format placeholders with
#' @return A formatted string
#'
#' @examples
#' \dontrun{
#' format_str("Format with {1} and {2}", "x", "y")
#' }
#'
#' @keywords internal
.format_str <- function(string, ...){
  args <- list(...)
  for (i in 1:length(args)){
    pattern <- paste("\\{", i, "}", sep="")
    replacement <- args[[i]]
    string <- gsub(pattern, replacement, string)
  }
  return(string)
}

#' @keywords internal
.check_hypeR_GEM_obj <- function(hypGEM_obj) {
  return(
    is(hypGEM_obj, "list") && all(
      purrr::map_vec(hypGEM_obj, \(sig) all(names(sig) == c("info", "data")))
    )
  )
}



#' Convert metabolite names to RefMet via Metabolomics Workbench
#'
#' Sends metabolite names to the Metabolomics Workbench RefMet converter and
#' returns standardized RefMet annotations and related metadata.
#'
#' @param DF A data.frame whose first column contains metabolite names.
#'
#' @return A data.frame with RefMet mapping results. Column names are
#'   normalized; includes `metabolite` and `refmet_name` when available.
#'
#' @export
#' @references Metabolomics Workbench RefMet:
#' \url{https://www.metabolomicsworkbench.org/databases/refmet/}
#'
#' @importFrom curl new_handle handle_setform handle_setopt curl_fetch_memory
#' @importFrom dplyr rename_with rename mutate
#' @importFrom stringr str_replace
refmet_convert <- function(DF) {
  # --- input checks ---
  if (!is.data.frame(DF) || ncol(DF) < 1L) {
    stop("`DF` must be a data.frame with at least one column.")
  }

  # Join names as newline-separated text (ensure character)
  mets <- paste(as.character(DF[[1L]]), collapse = "\n")

  # HTTP request
  h <- curl::new_handle()
  curl::handle_setform(h, metabolite_name = mets)
  curl::handle_setopt(h, timeout = 60)

  resp <- curl::curl_fetch_memory(
    "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_min.php",
    handle = h
  )
  if (is.null(resp$status_code) || resp$status_code != 200L) {
    stop("RefMet request failed with status code: ", resp$status_code %||% "unknown")
  }

  txt <- rawToChar(resp$content)
  if (!nzchar(txt)) return(data.frame())

  # Parse TSV robustly
  con <- textConnection(txt); on.exit(close(con), add = TRUE)
  refmet <- utils::read.delim(
    con, header = TRUE, sep = "\t", quote = "",
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # Drop all-blank rows
  if (nrow(refmet) == 0L) return(refmet)
  is_all_blank <- apply(refmet, 1L, function(r) all(is.na(r) | r == ""))
  refmet <- refmet[!is_all_blank, , drop = FALSE]

  # Normalize column names
  refmet <- dplyr::rename_with(refmet, ~ trimws(.x))|>
    dplyr::rename_with(refmet, ~ tolower(stringr::str_replace(.x, "[ ]+", "_")))

  # Safe renames if present
  if ("input_name" %in% names(refmet)) {
    refmet <- dplyr::rename(refmet, metabolite = "input_name")
  }
  if ("standardized_name" %in% names(refmet)) {
    refmet <- dplyr::rename(refmet, refmet_name = "standardized_name")
  }
  if ("refmet_name" %in% names(refmet)) {
    refmet <- dplyr::mutate(refmet, refmet_name = ifelse(refmet_name == "-", NA_character_, refmet_name))
  }

  return(refmet)
}

# null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b
