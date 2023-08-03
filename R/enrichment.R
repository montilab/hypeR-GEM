#' Title Enrichment analysis
#'
#' @param hypeR_GEM_obj an hypeR_GEM object, the output of "signature2gene" function
#' @param genesets a list of genesets
#' @param method enrichment method
#' @param weights the column name in the gene_table of hypeR_GEM_obj that represents the weight of each gene
#' @param background background parameter of hypergeometric test

#' @importFrom stats phyper p.adjust

#' @return a list of type of statistical test and data
#' @export
enrichment <- function(hypeR_GEM_obj,
                       genesets,
                       genesets_name = 'unknown',
                       method=c('unweighted','weighted','kstest'),
                       weights = 'one_minus_p_value',
                       background=1234567){

  # Default arguments
  method <- match.arg(method)

  ## multi_hyp_GEM
  if(method == "unweighted"){
    if(length(hypeR_GEM_obj$gene_tables)==1){
      hyp_GEM_enrichment <- .hyper_enrichment(hypeR_GEM_obj$gene_tables[[1]]$symbol,
                                              genesets = genesets,
                                              genesets_name = genesets_name,
                                              background = background)
    }else{
      hyp_GEM_enrichment <- lapply(hypeR_GEM_obj$gene_tables,
                                   function(x){return (x$symbol)}) %>%
        lapply(., .hyper_enrichment, genesets = genesets, genesets_name = genesets_name, background = background)
    }
  }



  ## multi_hyp_GEM
  if(method == "weighted"){
    if(length(hypeR_GEM_obj$gene_tables)==1){
      hyp_GEM_enrichment <- hypeR_GEM_obj$gene_tables[[1]] %>%
        dplyr::select(symbol,!!as.name(weights)) %>%
        tibble::deframe() %>%
        .weighted_hyper_enrichment(., genesets = genesets, genesets_name = genesets_name, background = background)
    }else{
      hyp_GEM_enrichment <- lapply(hypeR_GEM_obj$gene_tables,
                                   function(x){x <- x %>%
                                     dplyr::select(symbol,!!as.name(weights)) %>%
                                     tibble::deframe()}) %>%
        lapply(., .weighted_hyper_enrichment, genesets = genesets, genesets_name = genesets_name, background = background)
    }
  }

  return(hyp_GEM_enrichment)
}


#' Title Hypergeometric test
#'
#' @param unweighted_signature  a character vector of signature
#' @param genesets a list of genesets
#' @param background background parameter of hypergeometric test

#' @importFrom stats phyper p.adjust
#' @importFrom magrittr %>%

#' @return a list
#' @keywords internal
.hyper_enrichment <- function(unweighted_signature,
                              genesets,
                              genesets_name='unknown',
                              background=1234567){

  if (!is(genesets, "list")) stop("Expected genesets to be a list of genesets\n")

  signature <- unique(unweighted_signature)
  genesets <- lapply(genesets, unique)

  ## parameters
  signature_found <- signature[signature %in% unique(unlist(genesets))]
  hits <- sapply(genesets, function(x, y) length(intersect(x, y)), signature_found)
  drawn <- length(signature)
  n_genesets <- sapply(genesets, length)
  left <- background-n_genesets

  # Hypergeometric test
  pvals <- suppressWarnings(stats::phyper(q=hits-1,
                                          m=n_genesets,
                                          n=left,
                                          k=drawn,
                                          lower.tail=FALSE))

  # Format data
  data <- data.frame(label=names(genesets),
                     pval=signif(pvals, 2),
                     fdr=signif(stats::p.adjust(pvals, method="fdr"), 2),
                     signature=length(signature),
                     geneset=n_genesets,
                     overlap=hits,
                     weighted_overlap = hits,
                     background=background,
                     hits=sapply(genesets, function(x, y) paste(intersect(x, y), collapse=';'), signature_found),
                     stringsAsFactors=FALSE)

  return(list(info=list(Test = "Hypergeometric test",
                        Signature_size = length(signature),
                        Genesets = genesets_name,
                        Background = background),
              data=data))
}


#' Title Weighted hypergeometric test
#'
#' @param weighted_signature named vector, name = gene symbol, value = weight and should be between 0 and 1
#' @param genesets a list of genesets
#' @param background background parameter of hypergeometric test

#' @importFrom stats phyper p.adjust
#' @importFrom magrittr %>%

#' @return a list
#' @keywords internal
.weighted_hyper_enrichment <- function(weighted_signature,
                                      genesets,
                                      genesets_name='unknown',
                                      background=1234567){

  if (!is(weighted_signature, "vector")) stop("Expected signature to be a vector of symbols\n")
  if (!is(genesets, "list")) stop("Expected genesets to be a list of genesets\n")

  if(max(weighted_signature) > 1 | min(weighted_signature) < 0) stop("All weights should be between 0 and 1\n")

  if(is.null(names(weighted_signature))) stop("Expected signature to be a named vector\n")

  ## parameters
  signature <- unique(names(weighted_signature))
  genesets <- lapply(genesets,unique)
  signature_found <- signature[signature %in% unique(unlist(genesets))]


  ## weighted hits(overlaps)
  weighted_hits <- lapply(genesets, function(x, y) intersect(x, y), signature_found) %>%
    lapply(., function(x,y) y[x], weighted_signature) %>%
    lapply(., sum) %>%
    unlist(.) %>%
    round(.,0)

  ## For each geneset, all genes have a weight = 1
  n_genesets <- lengths(genesets)
  weighted_drawn <- round(sum(weighted_signature), 0)
  weighted_left <- background-n_genesets

  ## Hypergeometric test
  pvals <- suppressWarnings(stats::phyper(q=weighted_hits-1,
                                          m=n_genesets,
                                          n=weighted_left,
                                          k=weighted_drawn,
                                          lower.tail=FALSE))
  # Format data
  data <- data.frame(label=names(genesets),
                     pval=signif(pvals, 2),
                     fdr=signif(stats::p.adjust(pvals, method="fdr"), 2),
                     signature=length(signature),
                     geneset=n_genesets,
                     overlap=sapply(genesets, function(x, y) length(intersect(x, y)), signature_found),
                     weighted_overlap=weighted_hits,
                     background=background,
                     hits=sapply(genesets, function(x, y) paste(intersect(x, y), collapse=';'), signature_found),
                     stringsAsFactors=FALSE)

  return(list(info=list(Test = "Weighted hypergeometric test",
                        Signature_size = length(signature),
                        Genesets = genesets_name,
                        Background = background),
              data = data))
}
