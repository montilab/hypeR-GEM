#' Enrichment analysis
#'
#' @param hypeR_GEM_obj an hypeR_GEM object, the output of "signature2gene" function
#' @param genesets a list of genesets
#' @param genesets_name name of the geneset,e.g "KEGG"
#' @param method enrichment method
#' @param weights the column name in the gene_table of hypeR_GEM_obj that represents the weight of each gene
#' @param min_metabolite minimum number of metabolite that drives the enrichment of a pathway
#' @param background background parameter of hypergeometric test


#' @return a list of type of statistical test and data

#' @importFrom rlang .data
#' @import methods utils
#' @export
enrichment <- function(hypeR_GEM_obj,
                       genesets,
                       genesets_name = 'unknown',
                       method=c('unweighted','weighted'),
                       weights = 'one_minus_fdr',
                       min_metabolite = 0,
                       background=1234567){

  ## check
  if(!is.list(hypeR_GEM_obj)) stop("hypeR_GEM_obj must be a list object!\n")
  if(!is.list(genesets)) stop("genesets must be a list object!\n")
  if(!(any(names(hypeR_GEM_obj) %in% c("mapped_metabolite_signatures", "gene_tables")))) stop("element names in 'hypeR_GEM_obj' must be contain 'mapped_metabolite_signatures' and 'gene_tables'")
  if(min_metabolite < 0) stop("'min_metabolite' must be non-negative integer")

  weights_in_mapped_gene_tables <- all(lapply(hypeR_GEM_obj$gene_tables,colnames) %>%
        lapply(., function(x){return(weights %in% x)}) %>%
        unlist())

  if(!all(weights_in_mapped_gene_tables)) stop("'weights' must be one of the column in all element of hypeR_GEM_obj$gene_tables")

  # Default arguments
  method <- match.arg(method)
  min_metabolite <- floor(min_metabolite)


  ## unweighted hypergeometric
  if(method == "unweighted"){
    ## single or multiple signatures
    if(length(hypeR_GEM_obj$gene_tables)==1){
      hypeR_GEM_enrichments <- .hyper_enrichment(hypeR_GEM_obj$gene_tables[[1]],
                                                genesets = genesets,
                                                genesets_name = genesets_name,
                                                min_metabolite =  min_metabolite,
                                                background = background)
    }else{
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
    ## single or multiple signatures
    if(length(hypeR_GEM_obj$gene_tables)==1){
      hypeR_GEM_enrichments <- .weighted_hyper_enrichment(hypeR_GEM_obj$gene_tables[[1]],
                                                       genesets = genesets,
                                                       genesets_name = genesets_name,
                                                       min_metabolite= min_metabolite,
                                                       weights = weights,
                                                       background = background)
    }else{
      hypeR_GEM_enrichments <- lapply(hypeR_GEM_obj$gene_tables,
                                   .weighted_hyper_enrichment,
                                   genesets = genesets,
                                   genesets_name = genesets_name,
                                   min_metabolite = min_metabolite,
                                   weights = weights,
                                   background = background)
    }
  }

  return(hypeR_GEM_enrichments)
}


#' Title Hypergeometric test
#'
#' @param gene_table the gene table from hypeR-GEM object
#' @param genesets a list of genesets
#' @param genesets_name name of the geneset,e.g "KEGG"
#' @param min_metabolite minimum number/ratio of metabolite that drives the enrichment
#' @param background background parameter of hypergeometric test

#' @import methods utils
#' @importFrom stats phyper p.adjust
#' @importFrom magrittr %>%
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
#' @param background background parameter of hypergeometric test

#' @import methods utils
#' @importFrom stats phyper p.adjust
#' @importFrom magrittr %>%
#' @importFrom stringr str_count
#' @importFrom dplyr filter mutate arrange
#' @importFrom rlang .data

#' @return a list
#' @keywords internal
.weighted_hyper_enrichment <- function(gene_table,
                                       genesets,
                                       genesets_name='unknown',
                                       min_metabolite,
                                       weights = 'one_minus_fdr',
                                       background=1234567){


  if (!is(genesets, "list")) stop("Expected genesets to be a list of genesets\n")
  if(!(weights %in% colnames(gene_table))) stop("Gene weights must be specified in a colnmae of gene_table\n")

  weighted_signature <- gene_table %>%
    dplyr::select(symbol,!!as.name(weights)) %>%
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
