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
