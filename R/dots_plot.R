#' Title Create a dot plot
#'
#' @param hyp_enrichment_data
#' @param top
#' @param abrv
#' @param size_by
#' @param pval_cutoff
#' @param fdr_cutoff
#' @param val
#' @param title
#'
#' @return A ggplot object
#' @export
dots_plot <- function(hyp_enrichment_data,
                       top=20,
                       abrv=50,
                       size_by=c("genesets", "significance", "none"),
                       pval_cutoff=1,
                       fdr_cutoff=1,
                       val=c("fdr", "pval"),
                       title="") {

  # Default arguments
  val <- match.arg(val)
  size_by <- match.arg(size_by)

  # Subset results
  df <- hyp_enrichment_data %>%
    dplyr::filter(pval <= pval_cutoff) %>%
    dplyr::filter(fdr <= fdr_cutoff) %>%
    purrr::when(!is.null(top) ~ head(., top), ~ .)

  # Handle empty dataframes
  if (nrow(df) == 0) return(ggempty())

  # Plotting variables
  df$significance <- df[,val]
  df$size <- 1

  if (size_by == "significance") {
    df$size <- df$significance
  }
  if (size_by == "genesets") {
    df$size <- df$geneset
  }

  # Order by significance value
  df <- df[order(-df[,val]),]

  # Abbreviate labels
  label.abrv <- substr(df$label, 1, abrv)
  if (any(duplicated(label.abrv))) {
    stop("Non-unique labels after abbreviating")
  } else {
    df$label.abrv <- factor(label.abrv, levels=label.abrv)
  }

  if (val == "pval") {
    color.label <- "P-Value"
  }
  if (val == "fdr") {
    color.label <- "FDR"
  }

  p <- ggplot(df, aes(x=label.abrv, y=significance, color=significance, size=size)) +
    geom_point() +
    labs(title=title, y=color.label, color=color.label) +
    scale_color_continuous(low="#E53935", high="#114357", guide=guide_colorbar(reverse=TRUE)) +
    coord_flip() +
    scale_y_continuous(trans=ggforce::trans_reverser("log10")) +
    geom_hline(yintercept=0.05, linetype="dotted") +
    theme(plot.title=element_text(hjust=0.5),
          axis.title.y=element_blank())

  if (size_by == "none") {
    p <- p + guides(size="none")
  }
  if (size_by == "significance") {
    p <- p + scale_size_continuous(trans=ggforce::trans_reverser("log10")) + labs(size="Significance")
  }
  if (size_by == "genesets") {
    p <- p + scale_size_continuous(trans=scales::log10_trans()) + labs(size="Genesets\nSize")
  }

  return(p)


}

