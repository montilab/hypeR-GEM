#' Title Map metabolic signatures to enzyme-coding genes
#'
#' @param signatures a list of metabolic signatures, each element is a data frame which has to contain a column with the same name as "reference_key"
#' @param species the species of GEM model, if species == 'Other', a user-define GEM_table must be provided
#' @param directed logical parameter, if TRUE, map metabolites to reactions where these metabolites are product only
#' @param tables species == 'Other', a user-define GEM_table must be provided
#' @param merge merged metabolites in different compartment, for example, MAM001_c, MAM001_e, ... -> MAM001
#' @param reference_key the column name in the reference table which represents the standardized names (e.g. "refmet_name")
#' @param ensemble_id if the genes in GEM_tables$gene_df is given by Ensemble, then ensemble_id = TRUE, otherwise, ensemble_id = FALSE
#' @param promiscuous_threshold gene association threshold of promiscuous metabolite
#' @param background the background of the gene-specific hypergeometric test, default = NULL = number of merged metabolites
#' @return a list containing two element: "mapped_metabolite_signatures" and "gene_tables" for each signature

#' @import methods utils
#' @importFrom magrittr %>% is_greater_than
#' @importFrom dplyr filter select mutate left_join rename distinct pull
#' @importFrom Matrix Matrix
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom gprofiler2 gconvert

#' @export
signature2gene <- function(signatures,
                           species = c("Human", "Mouse", "Rat", "Zebrafish", "Worm", "Other"),
                           directed = FALSE,
                           tables = NULL,
                           merge = TRUE,
                           reference_key="refmet_name",
                           promiscuous_threshold = NULL,
                           ensemble_id = TRUE,
                           background = NULL){

  ## check both signatures and meta_df contain "reference_key"
  if(!is.list(signatures)) stop("signatures must be a list object!\n")
  if(is.null(names(signatures))) stop("Expected signature to be a named list\n")
  if(!all(sapply(signatures, function(x,reference_key){return(reference_key %in% colnames(x))}, reference_key=reference_key))) stop("reference key does not in signatures' column names!\n")

  # Default arguments
  species <- match.arg(species)

  # load tables
  if(species == "Human" & !directed){GEM_tables <- hypeR.GEM::Human_GEM_tables}
  if(species == "Human" &  directed){GEM_tables <- hypeR.GEM::Human_GEM_tables_directed}

  if(species == "Mouse" & !directed){GEM_tables <- hypeR.GEM::Mouse_GEM_tables}
  if(species == "Mouse" &  directed){GEM_tables <- hypeR.GEM::Mouse_GEM_tables_directed}

  if(species == "Rat" & !directed){GEM_tables <- hypeR.GEM::Rat_GEM_tables}
  if(species == "Rat" &  directed){GEM_tables <- hypeR.GEM::Rat_GEM_tables_directed}


  if(species == "Zebrafish" & !directed){GEM_tables <- hypeR.GEM::Zebrafish_GEM_tables}
  if(species == "Zebrafish" & directed){GEM_tables <- hypeR.GEM::Zebrafish_GEM_tables_directed}

  if(species == "Worm" & !directed){GEM_tables <- hypeR.GEM::Worm_GEM_tables}
  if(species == "Worm" & directed){GEM_tables <- hypeR.GEM::Worm_GEM_tables_directed}


  ## user-defined GEM table
  if(species == "Other" & is.null(tables)) stop("GEM tables must be provided as a list for an unspecified species!")

  if(species == "Other" & !is.null(tables)){
    if(!(reference_key %in% colnames(GEM_tables[['meta_df']]))) stop("reference key does not in GEM's metabolite dataframe column names!\n")
    GEM_tables <- tables
  }

  ## if compartments are merged
  if(merge){
    meta_df <- GEM_tables$meta_df_merged
    m2g <- GEM_tables$m2g_merged
    g2m <- GEM_tables$g2m_merged
    g2r <- GEM_tables$g2r
  }else{
    meta_df <- GEM_tables$meta_df
    m2g <- GEM_tables$m2g
    g2m <- GEM_tables$g2m
    g2r <- GEM_tables$g2r
  }



  ## map metabolic signatures to metabolites in the GEM via "reference_key"
  mapped_metabolite_signatures <- lapply(signatures, .extract_standardized_name, reference_key) %>%
    lapply(., .extract_sig_meta, meta_df, reference_key) %>%
    lapply(., .gene_association, m2g, reference_key) %>%
    lapply(., function(x){x %>%
        dplyr::distinct(!!as.name(reference_key), .keep_all = TRUE)})

  ## remove promiscuous metabolites
  if(!is.null(promiscuous_threshold)){
    promiscuous_met <- lapply(mapped_metabolite_signatures, function(x){
      x %>%
        dplyr::filter(gene_association > promiscuous_threshold ) %>%
        dplyr::select(!!as.name(reference_key), gene_association)
    }) %>%
      do.call(rbind, .) %>%
      magrittr::set_rownames(., NULL) %>%
      dplyr::distinct()

    mapped_metabolite_signatures <- lapply(mapped_metabolite_signatures,function(x){
      x %>%
        dplyr::filter(!(!!as.name(reference_key) %in%  promiscuous_met[[reference_key]]))
    })
  }



  ## map metabolic signatures to ranked genes
  ### default background = # of unique metabolites in GEM
  if(is.null(background)){
    background <- GEM_tables[['meta_df_merged']] %>%
      nrow(.)
  }

  gene_tables <- lapply(mapped_metabolite_signatures,
                        .meta2gene,
                        meta_df = meta_df,
                        m2g=m2g,
                        g2r=g2r,
                        g2m=g2m,
                        background=background,
                        reference_key=reference_key)

  if(ensemble_id){
    ## map from ENSEMBL id to gene symbol
    gene_map <- gprofiler2::gconvert(query = unlist(unique(lapply(gene_tables,rownames))),
                                     organism = "hsapiens",
                                     target = "ENSG",
                                     filter_na = TRUE)

    ## map ENSEMBL id in gene_tables to gene symbol
    gene_tables <- lapply(gene_tables, function(x){x <- x %>%
      dplyr::filter(row.names(.) %in% gene_map[['input']]) %>%
      tibble::rownames_to_column(var='input') %>%
      dplyr::left_join(gene_map, by='input') %>%
      dplyr::distinct(input,.keep_all = TRUE) %>%
      dplyr::rename(symbol = name,
                    name = input) %>%
      dplyr::select(name,
                    symbol,
                    associated_reactions,
                    total_association,
                    signature_association,
                    signature_size,
                    background,
                    p_value,
                    fdr,
                    one_minus_p_value,
                    one_minus_fdr,
                    associated_metabolites)})
  }else{
    gene_tables <- lapply(gene_tables, function(x){x %>% tibble::rownames_to_column(var='symbol')})
  }

  if(!is.null(promiscuous_threshold)){
    return(list(promiscuous_met = promiscuous_met,
                mapped_metabolite_signatures = mapped_metabolite_signatures,
                gene_tables = gene_tables))
  }else{
    return(list(mapped_metabolite_signatures = mapped_metabolite_signatures,
                gene_tables = gene_tables))
  }
}


#' Title Extract the standardized names from the signature data frames
#'
#' @param signature_df signature data frame
#' @param reference_key the key which is used to map the signature and metabolites in the GEM

#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull

#' @return a character vector contains standardized names of signatures
#' @keywords internal
.extract_standardized_name <- function(signature_df, reference_key="refmet_name"){

  names <- signature_df %>%
    dplyr::filter(!is.na(!!as.name(reference_key))) %>%
    dplyr::pull(!!as.name(reference_key))

  return(names)
}


#' Title Subset the GEM metabolite data frame
#'
#' @param signature a character vector contains standardized names of signatures
#' @param GEM_meta_df the GEM metabolite data frame
#' @param reference_key the key which is used to map the signature and metabolites in the GEM

#' @importFrom magrittr %>%
#' @importFrom dplyr filter

#' @return a subset of the GEM metabolite data frame
#' @keywords internal
.extract_sig_meta <- function(signature,
                              GEM_meta_df,
                              reference_key="refmet_name"){

  df <-  GEM_meta_df %>%
    dplyr::filter(!!as.name(reference_key) %in%  signature)
  return(df)
}

#' Title Gene association of each signature
#'
#' @param mapped_metabolite_signature one element in "mapped_metabolite_signatures"
#' @param m2g association list
#' @param reference_key the key which is used to map the signature and metabolites in the GEM

#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select

#' @return a updated 'mapped_metabolite_signature' with "gene_association" column
#' @keywords internal
.gene_association <- function(mapped_metabolite_signature,
                             m2g,
                             reference_key = "refmet_name"){

  mapped_metabolite_signature <- mapped_metabolite_signature %>%
    dplyr::mutate(gene_association  = lengths(m2g[mapped_metabolite_signature$name])) %>%
    dplyr::select(name, fullname, !!as.name(reference_key),  gene_association, everything())
}

#' Title Map metabolic signature to ranked genes
#'
#' @param signature metabolite signature
#' @param m2g association list
#' @param g2r association list
#' @param g2m association list
#' @param background background for gene-specific hypergeometric test
#' @param reference_key the key which is used to map the signature and metabolites in the GEM

#' @importFrom magrittr %>%
#' @importFrom dplyr filter select pull
#' @importFrom stats phyper p.adjust

#' @return a data frame
#' @keywords internal
.meta2gene <- function(signature,
                       meta_df,
                       m2g,
                       g2r,
                       g2m,
                       background,
                       reference_key='refmet_name'){

  ## map metabolite signatures to genes
  mapped_genes <- m2g[signature$name] %>%
    unlist(.) %>%
    unique(.)

  ## compute statistics
  associated_reactions <- lengths(g2r[mapped_genes])

  total_association <- lengths(g2m[mapped_genes])

  signature_association <- lengths(lapply(g2m[mapped_genes], function(x,y){return(intersect(x,y))}, signature$name))

  ## organize statistics in a data frame
  gene_df <- data.frame(associated_reactions = associated_reactions,
                        total_association = total_association,
                        signature_association = signature_association,
                        signature_size = rep(nrow(signature), length(mapped_genes)),
                        background = rep(background, length(mapped_genes))) %>%
    dplyr::mutate(p_value = suppressWarnings(stats::phyper(q=signature_association,
                                                           m=signature_size,
                                                           n=background - signature_size,
                                                           k=total_association,
                                                           lower.tail=FALSE)),
                  fdr = stats::p.adjust(p_value, method="fdr"),
                  one_minus_p_value = 1-p_value,
                  one_minus_fdr = 1-fdr)

  ## associated metabolites in signature of each gene
  gene_df$associated_metabolites <- lapply(g2m[mapped_genes], function(x,y){return(intersect(x,y))}, signature$name) %>%
    lapply(., function(x){meta_df %>%
        dplyr::filter(name %in% x) %>%
        dplyr::pull(!!as.name(reference_key))}) %>%
    lapply(., function(x){return(paste(x,collapse = ";"))}) %>%
    unlist(.)

  return(gene_df)
}
