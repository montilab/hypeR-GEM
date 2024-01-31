#' Title Reactable table for hypeR_GEM_enrichment object
#'
#' @param hypeR_GEM_obj A list of enrichment results from single/multiple signatures, output of "hypeR.GEM::enrichment()"
#' @param fdr_cutoff fdr threshold for geneset enrichment

#' @importFrom reactable reactable colDef
#' @importFrom htmltools div

#' @return a nested reactable
#' @export
rctbls <- function(
    hypeR_GEM_obj,
    fdr_cutoff = 0.05
) {
  ## input checks
  stopifnot( .check_hypeR_GEM_obj(hypeR_GEM_obj) )

  ## filter out non-enriched geneset
  hypeR_GEM_obj <- hypeR_GEM_obj |>
    purrr::map(\(ls) {
      ls$data <- ls$data |> dplyr::filter(fdr < fdr_cutoff)
      return(ls)
    })
  ## create outer reactable
  outer_df <- data.frame(
    signature = names(hypeR_GEM_obj),
    size = purrr::map_vec(
      hypeR_GEM_obj, \(x) x$info$Signature_size, .ptype = integer(1)
    ),
    enriched = purrr::map_vec(
      hypeR_GEM_obj, \(x) nrow(x$data), .ptype = integer(1)
    ),
    gsets = purrr::map_vec(
      hypeR_GEM_obj, \(x) x$info$Genesets, .ptype = character(1)
    ),
    bg = purrr::map_vec(
      hypeR_GEM_obj, \(x) x$info$Background, .ptype = integer(1)
    )
  )
  ## inner reactable
  tbl <- reactable(outer_df,
    showPageSizeOptions = FALSE,
    onClick = "expand",
    resizable = TRUE,
    rownames = FALSE,
    defaultColDef = colDef(headerClass = "rctbl-outer-header"),
    columns = list(
      signature = colDef(name = "Signature", minWidth = 300),
      size = colDef(name = "Signature Size"),
      enriched = colDef(name = "Enriched Genesets"),
      gsets = colDef(name = "Genesets"),
      bg = colDef(name = "Background")
    ),
    details = function(index) {
      df <- .rctbl(hypeR_GEM_obj[[index]], type = "inner", fdr_cutoff = fdr_cutoff)
    },
    wrap = FALSE,
    class = "rctbl-outer-tbl",
    rowStyle = list(cursor = "pointer")
  )
  dat <- htmltools::div(class = "rctbl-outer-obj", tbl)
  return(dat)
}
#' Title Reactable table for hypeR_GEM_enrichment object
#'
#' @param hypeR_GEM_enrichment A hypeR_GEM_enrichment from a single signature
#' @param type Use style class for outer or inner tables

#' @importFrom reactable reactable colDef
#' @importFrom htmltools div tagAppendChild
#' @importFrom dplyr select

#' @return a reactable
#' @keywords internal
.rctbl <- function(
    hypeR_GEM_enrichment,
    type = c("inner", "outer"),
    fdr_cutoff = 0.05
)
{
  type <- match.arg(type)
  class.obj <- .format_str("rctbl-{1}-obj", type)
  class.tbl <- .format_str("rctbl-{1}-tbl", type)
  class.header <- .format_str("rctbl-{1}-header", type)

  df <- hypeR_GEM_enrichment$data |>
    dplyr::select(c("label", "pval", "fdr", "geneset", "overlap",
                    "weighted_overlap", "gene_hits", "metabolite_hits"))

  tbl <- reactable(df,
    rownames = FALSE,
    resizable = TRUE,
    searchable = TRUE,
    showPageSizeOptions = FALSE,
    compact = TRUE,
    defaultColDef = colDef(headerClass = class.header),
    columns = list(
      label = colDef(name = "Label", minWidth = 300),
      pval = colDef(name = "P-Value"),
      fdr = colDef(name = "FDR"),
      geneset = colDef(name = "Geneset Size"),
      overlap = colDef(name = "Overlap"),
      weighted_overlap = colDef(name = "Weighted overlap"),
      gene_hits = colDef(name = "Gene Hits"),
      metabolite_hist = colDef(name = "Metabolite Hits")
    ),
    class = class.tbl
  )
  dat <- htmltools::div(class = class.obj, tbl)

  return(dat)
}



