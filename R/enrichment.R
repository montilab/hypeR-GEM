#' Enrichment analysis
#'
#' @param hypeR_GEM_obj an hypeR_GEM object, the output of "signature2gene" function
#' @param genesets a list of genesets
#' @param genesets_name name of the geneset,e.g "KEGG"
#' @param method enrichment method
#' @param weighted_by the column name in the gene_table of hypeR_GEM_obj that represents the weight of each gene
#' @param a proportional to the smoothness of the sigmoid function, default a = 1
#' @param b the half-point threshold of the sigmoid function, default b = -1 (p-value = 0.1 as half-point)
#' @param sigmoid_transformation logical; when method == "weighted", if TRUE apply .sigmoid_transformation() to `weighted_by`; if FALSE use the raw values
#' @param min_metabolite minimum number of metabolite that drives the enrichment of a pathway
#' @param background background parameter of hypergeometric test
#'
#' @return a list of type of statistical test and data
#' @importFrom rlang .data
#' @import methods utils
#' @export
enrichment <- function(hypeR_GEM_obj,
                       genesets,
                       genesets_name = 'unknown',
                       method = c('unweighted','weighted'),
                       weighted_by = 'one_minus_fdr',
                       a =  1,
                       b = -1,
                       sigmoid_transformation = TRUE,
                       min_metabolite = 0,
                       background = 1234567){

  ## check
  if(!is.list(hypeR_GEM_obj)) stop("hypeR_GEM_obj must be a list object!\n")
  if(!is.list(genesets)) stop("genesets must be a list object!\n")
  if(!all(c("mapped_metabolite_signatures","gene_tables") %in% names(hypeR_GEM_obj)))
    stop("element names in 'hypeR_GEM_obj' must contain 'mapped_metabolite_signatures' and 'gene_tables'")
  if(min_metabolite < 0) stop("'min_metabolite' must be non-negative integer")

  # Default arguments
  method <- match.arg(method)
  min_metabolite <- floor(min_metabolite)

  ## only validate weight column if weighted method is chosen
  if (method == "weighted") {
    col_ok <- vapply(hypeR_GEM_obj$gene_tables,
                     function(tbl) weighted_by %in% colnames(tbl),
                     logical(1))
    if (!all(col_ok)) {
      stop(sprintf("'%s' must be a column in every element of hypeR_GEM_obj$gene_tables", weighted_by))
    }
  }

  ## unweighted hypergeometric
  if(method == "unweighted"){
    if(length(hypeR_GEM_obj$gene_tables)==1){
      hypeR_GEM_enrichments <- .hyper_enrichment(hypeR_GEM_obj$gene_tables[[1]],
                                                 genesets = genesets,
                                                 genesets_name = genesets_name,
                                                 min_metabolite =  min_metabolite,
                                                 background = background)
    } else {
      hypeR_GEM_enrichments <- lapply(hypeR_GEM_obj$gene_tables,
                                      .hyper_enrichment,
                                      genesets = genesets,
                                      genesets_name = genesets_name,
                                      min_metabolite =  min_metabolite,
                                      background = background)
    }
  }

  ## weighted hypergeometric
  if(method == "weighted"){
    if(length(hypeR_GEM_obj$gene_tables)==1){
      hypeR_GEM_enrichments <- .weighted_hyper_enrichment(hypeR_GEM_obj$gene_tables[[1]],
                                                          genesets = genesets,
                                                          genesets_name = genesets_name,
                                                          min_metabolite = min_metabolite,
                                                          weighted_by = weighted_by,
                                                          a = a,
                                                          b = b,
                                                          sigmoid_transformation = sigmoid_transformation,
                                                          background = background)
    } else {
      hypeR_GEM_enrichments <- lapply(hypeR_GEM_obj$gene_tables,
                                      .weighted_hyper_enrichment,
                                      genesets = genesets,
                                      genesets_name = genesets_name,
                                      min_metabolite = min_metabolite,
                                      weighted_by = weighted_by,
                                      a = a,
                                      b = b,
                                      sigmoid_transformation = sigmoid_transformation,
                                      background = background)
    }
  }

  return(hypeR_GEM_enrichments)
}



#' Hypergeometric test
#'
#' @param gene_table the gene table from hypeR-GEM object
#' @param genesets a list of genesets
#' @param genesets_name name of the geneset,e.g "KEGG"
#' @param min_metabolite minimum number/ratio of metabolite that drives the enrichment
#' @param background background parameter of hypergeometric test

#' @import methods utils
#' @importFrom stats phyper p.adjust
#' @importFrom magrittr %>%
#' @importFrom tibble deframe
#' @importFrom stringr str_count
#' @importFrom dplyr filter mutate arrange
#' @importFrom rlang .data

#' @return a list
#' @keywords internal
.hyper_enrichment <- function(gene_table,
                              genesets,
                              genesets_name='unknown',
                              min_metabolite,
                              background=1234567){

  if (!is(genesets, "list")) stop("Expected genesets to be a list of genesets\n")

  signature <- unique(gene_table$symbol)
  genesets <- lapply(genesets, unique)

  ## parameters
  signature_found <- signature[signature %in% unique(unlist(genesets))]
  hits <- sapply(genesets, function(x, y) length(intersect(x, y)), signature_found)
  drawn <- length(signature)
  n_genesets <- sapply(genesets, length)
  left <- background-n_genesets

  ## record metabolite signature size
  metabolite_signature_size <- unique(gene_table$signature_size)


  ## Associated metabolites
  hitted_genes <- lapply(genesets, function(x, y){intersect(x, y)}, signature_found)
  metabolite_hits <- lapply(hitted_genes, function(x){gene_table %>%
      dplyr::filter(symbol %in% x) %>%
      dplyr::pull(associated_metabolites)}) %>%
    lapply(., strsplit, ";") %>%
    lapply(., unlist) %>%
    lapply(., unique) %>%
    lapply(., paste, collapse=";") %>%
    unlist()

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
                     background=background,
                     gene_hits=unlist(lapply(hitted_genes, paste, collapse=";")),
                     metabolite_hits = metabolite_hits,
                     stringsAsFactors=FALSE) %>%
    dplyr::mutate(num_met_hits = stringr::str_count(metabolite_hits, ";") + 1,
                  ratio_met_hits = round(num_met_hits/metabolite_signature_size, 3)) %>%
    dplyr::filter(num_met_hits >= min_metabolite) %>%
    dplyr::arrange(pval)


  return(list(info=list(Test = "Hypergeometric test",
                        Signature_size = length(signature),
                        Genesets = genesets_name,
                        Background = background),
                        data=data))
}


#' Weighted hypergeometric test
#'
#' @param gene_table the gene table from hypeR-GEM object
#' @param genesets a list of genesets
#' @param genesets_name name of the geneset,e.g "KEGG"
#' @param min_metabolite minimum number/ratio of metabolite that drives the enrichment
#' @param weighted_by the value used for (optional) sigmoid transformation
#' @param a proportional to the smoothness of the sigmoid function, default a = -1
#' @param b the half-point threshold of the sigmoid function, default b = -1 (p-value = 0.1 as half-point)
#' @param sigmoid_transformation logical; if TRUE apply .sigmoid_transformation to `weighted_by`, else use raw values
#' @param background background parameter of hypergeometric test
#'
#' @import methods utils
#' @importFrom stats phyper p.adjust
#' @importFrom magrittr %>%
#' @importFrom tibble deframe
#' @importFrom stringr str_count
#' @importFrom dplyr filter mutate arrange
#' @importFrom rlang .data
#' @return a list
#' @keywords internal
.weighted_hyper_enrichment <- function(gene_table,
                                       genesets,
                                       genesets_name='unknown',
                                       min_metabolite,
                                       weighted_by = 'fdr',
                                       a = -1,
                                       b = -1,
                                       sigmoid_transformation = TRUE,
                                       background=1234567){

  if (!is(genesets, "list")) stop("Expected genesets to be a list of genesets\n")
  if(!(weighted_by %in% colnames(gene_table))) stop("Gene weights must be specified in a column of gene_table\n")

  ## decide weights
  if (isTRUE(sigmoid_transformation)) {
    gene_table$weights <- .sigmoid_transformation(
      p = gene_table[[weighted_by]],
      a = a,
      b = b
    )
  } else {
    gene_table$weights <- as.numeric(gene_table[[weighted_by]])
  }

  if (anyNA(gene_table$weights)) {
    stop(sprintf("'%s' contains NA values after transformation; please clean or impute.", weighted_by))
  }

  weighted_signature <- gene_table %>%
    dplyr::select(symbol, weights) %>%
    tibble::deframe()

  if(max(weighted_signature) > 1 | min(weighted_signature) < 0) stop("All weights should be between 0 and 1\n")

  ## parameters
  signature <- unique(names(weighted_signature))
  genesets <- lapply(genesets,unique)
  signature_found <- signature[signature %in% unique(unlist(genesets))]

  ## record metabolite signature size
  metabolite_signature_size <- unique(gene_table$signature_size)


  ## weighted hits(overlaps)
  weighted_hits <- lapply(genesets, function(x, y) intersect(x, y), signature_found) %>%
    lapply(., function(x,y) y[x], weighted_signature) %>%
    lapply(., sum) %>%
    unlist() %>%
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
  ## Associated metabolites
  hitted_genes <- lapply(genesets, function(x, y){intersect(x, y)}, signature_found)
  metabolite_hits <- lapply(hitted_genes, function(x){gene_table %>%
      dplyr::filter(symbol %in% x) %>%
      dplyr::pull(associated_metabolites)}) %>%
    lapply(., strsplit, ";") %>%
    lapply(., unlist) %>%
    lapply(., unique) %>%
    lapply(., paste, collapse=";") %>%
    unlist()

  # Format data
  data <- data.frame(label=names(genesets),
                     pval=signif(pvals, 2),
                     fdr=signif(stats::p.adjust(pvals, method="fdr"), 2),
                     signature=length(signature),
                     geneset=n_genesets,
                     weighted_overlap=weighted_hits,
                     background=background,
                     gene_hits=unlist(lapply(hitted_genes, paste, collapse=";")),
                     metabolite_hits = metabolite_hits,
                     stringsAsFactors=FALSE) %>%
    dplyr::mutate(num_met_hits = stringr::str_count(metabolite_hits, ";") + 1,
                  ratio_met_hits = round(num_met_hits/metabolite_signature_size, 3)) %>%
    dplyr::filter(num_met_hits >= min_metabolite) %>%
    dplyr::arrange(pval)


  return(list(info=list(Test = "Weighted hypergeometric test",
                        Signature_size = length(signature),
                        Genesets = genesets_name,
                        Background = background),
                        data = data))
}


#' @keywords internal
.sigmoid_transformation <- function(p,
                                    a = 1, ## smoothness
                                    b = -1,   # p = 0.1 as half-point
                                    eps = 1e-12) {

  # clip to avoid log10(0)
  p <- pmax(pmin(p, 1 - eps), eps)

  w <- 1 / (1 + exp(a * (log10(p) - b)))
  return(w)
}
