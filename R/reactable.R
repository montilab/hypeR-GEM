
#' Title Reactable table for hypeR_GEM_enrichment object
#'
#' @param hypeR_GEM_enrichments A hypeR_GEM_enrichments from multiple signatures

#' @importFrom reactable reactable colDef
#' @importFrom htmltools div

#' @return a nested reactable
#' @export
rctbls <- function(
    hypeR_GEM_enrichments,
    fdr_cutoff = 0.05) 
{
  hypeR_GEM_enrichments <- lapply(hypeR_GEM_enrichments, function(x) {
    x$data <- x$data |>
      dplyr::filter(fdr < fdr_cutoff)
    return(x)
  })
  outer_df <- data.frame(
    signature = names(hypeR_GEM_enrichments),
    size = sapply(hypeR_GEM_enrichments, function(x) {
      x$info[["Signature_size"]]
    }),
    enriched = sapply(hypeR_GEM_enrichments, function(x) {
      nrow(x$data)
    }),
    gsets = sapply(hypeR_GEM_enrichments, function(x) {
      x$info[["Genesets"]]
    }),
    bg = sapply(hypeR_GEM_enrichments, function(x) {
      x$info[["Background"]]
    })
  )
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
      df <- .rctbl(hypeR_GEM_enrichments[[index]], type = "inner", fdr_cutoff = fdr_cutoff)
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
#' @keyword internal
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



