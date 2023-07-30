#' Title Enrichment analysis visualization
#'
#' @param hypeR_GEM_enrichment an hypeR_GEM object, the output of "signature2gene" function
#' @param top
#' @param abrv
#' @param size_by
#' @param pval_cutoff
#' @param fdr_cutoff
#' @param val
#' @param title
#'
#' @return a ggplot object
#' @export
enrichment_plot <- function(hypeR_GEM_enrichment,
         top=20,
         abrv=50,
         size_by=c("genesets", "significance", "none"),
         pval_cutoff=1,
         fdr_cutoff=1,
         val=c("fdr", "pval"),
         title="") {

  # Default arguments
  size_by <- match.arg(size_by)
  val = match.arg(val)

  if(c('info') %in% names(hypeR_GEM_enrichment) & c('data') %in% names(hypeR_GEM_enrichment)){
    p <- .dots_plot(hypeR_GEM_enrichment$data,
                   top=40,
                   abrv=50,
                   size_by=size_by,
                   fdr_cutoff=max_fdr,
                   val=val)
  }else{
    p <- .multi_dots_plot(hypeR_GEM_enrichment,
                    top=40,
                    abrv=50,
                    size_by=size_by,
                    fdr_cutoff=max_fdr,
                    val=val)
  }

  return(p)
}


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
#' @keywords internal
.dots_plot <- function(hyp_enrichment_data,
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

#' Title Create a multi dots plot
#'
#' @param multihyp_data
#' @param top
#' @param abrv
#' @param size_by
#' @param pval_cutoff
#' @param fdr_cutoff
#' @param val
#' @param title
#'

#' @keywords internal
.multi_dots_plot <- function(multihyp_data,
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

  # Count significant genesets across signatures
  multihyp_dfs <- lapply(multihyp_data, function(hyp_obj) {
    hyp_obj$data %>%
      dplyr::filter(pval <= pval_cutoff) %>%
      dplyr::filter(fdr <= fdr_cutoff) %>%
      dplyr::select(label)
  })

  # Take top genesets
  labels <- names(sort(table(unlist(multihyp_dfs)), decreasing=TRUE))
  if (!is.null(top)) labels <- head(labels, top)
  # Handle empty dataframes
  if (length(labels) == 0) return(ggempty())

  # Create a multihyp dataframe
  dfs <- lapply(multihyp_data, function(hyp_obj) {
    hyp_df <- hyp_obj$data
    hyp_df[hyp_df$label %in% labels, c("label", val), drop=FALSE]
  })

  df <- suppressWarnings(Reduce(function(x, y) merge(x, y, by="label", all=TRUE), dfs))
  colnames(df) <- c("label", names(dfs))
  rownames(df) <- df$label
  df <- df[rev(labels), names(dfs)] %>%
    tibble::rownames_to_column(var='label')

  ## Add 'geneset' column to df
  df <- lapply(multihyp_data, function(hyp_obj) {
    hyp_obj$data[, c("label", "geneset")]
  }) %>%
    do.call(rbind, .) %>%
    dplyr::distinct(label, .keep_all=TRUE) %>%
    dplyr::filter(label %in% df$label) %>%
    dplyr::right_join(., df, by='label') %>%
    tibble::column_to_rownames(var='label')

  # Abbreviate labels
  label.abrv <- substr(rownames(df), 1, abrv)
  if (any(duplicated(label.abrv))) {
    stop("Non-unique labels after abbreviating")
  } else {
    rownames(df) <- factor(label.abrv, levels=label.abrv)
  }

  if (val == "pval") {
    cutoff <- pval_cutoff
    color.label <- "P-Value"
  }
  if (val == "fdr") {
    cutoff <- fdr_cutoff
    color.label <- "FDR"
  }

  ## separate "geneset" column
  geneset_size <- df %>%
    tibble::rownames_to_column(var='label') %>%
    dplyr::select(label, geneset)

  df.melted <- reshape2::melt(as.matrix(df %>% select(-c('geneset'))))
  colnames(df.melted) <- c("label", "signature", "significance")
  df.melted <- df.melted %>%
    dplyr::left_join(geneset_size, by='label')


  ## assign "size"
  if (size_by == "significance") {
    df.melted$size <- df.melted$significance
  }

  if (size_by == "genesets") {
    df.melted$size <- df.melted$geneset
  }

  if (size_by == "none") {
    df.melted$size <- 1
  }

  p <- df.melted %>%
    dplyr::filter(significance <= cutoff) %>%
    ggplot(aes(x=signature, y=label, color=significance, size=size)) +
    geom_point()+
    scale_color_continuous(low="#114357", high="#E53935", trans=ggforce::trans_reverser("log10")) +
    theme(plot.title=element_text(hjust=0.5),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1))

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


