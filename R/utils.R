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
.check_hypeR_GEM_obj <- function(hypeR_GEM_obj) {
  return(
    is(hypeR_GEM_obj, "list") && all(
      purrr::map_vec(hypeR_GEM_obj, \(sig) all(names(sig) == c("info", "data")))
    )
  )
}
