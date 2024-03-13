#' Title extract tables from a given genome scale metabolic model(GEM) of ".mat" format
#'
#' @param con a connection or the path of the GEM, which should be in .mat format
#' @param directional logical parameter, if TRUE, map metabolites to reactions where these metabolites are product only
#' @return a list of tables

#' @import methods utils
#' @importFrom R.matlab readMat
#' @importFrom magrittr %>%
#' @importFrom Matrix Matrix summary
#' @importFrom dplyr arrange filter select pull


#' @export
extract_tables <- function(con,
                           directional=FALSE){

  ## load GEM and extract data
  GEM <- R.matlab::readMat(con=file.path(con))
  data <- GEM[[names(GEM)]][,,1]

  ## create basic data frames of metabolites, reactions, and genes
  #compos <- sapply(data$comps[,1],unlist)
  #compartment <- compos[sapply(data$metComps[,1],unlist)]
  meta_df <- data.frame(name = sapply(data$mets[,1],unlist),
                        fullname = sapply(data$metNames[,1],unlist))

  reaction_df <- data.frame(name = sapply(data$rxns[,1],unlist))

  gene_df <- data.frame(name = unlist(sapply(data$genes,unlist)))

  ## map metabolites to reactions, row = reactions, column = metabolites
  r2m_matrix <- matrix(0,nrow=nrow(reaction_df),
                       ncol=nrow(meta_df),
                       dimnames = list(reaction_df$name,meta_df$name))

  ### pull out the S matrix, Matrix::summary() converts S@i, S@j from 0-based to 1-based
  S <- as(as(as(data$S, "dMatrix"), "generalMatrix"), "TsparseMatrix") %>%
    Matrix::summary(.) %>%
    as.data.frame(.) %>%
    dplyr::arrange(i)

  if(directional){
    ## Metabolite X maps to reactions where X is the product
    for(m in 1:ncol(r2m_matrix)){
      reactions <- S %>%
        dplyr::filter(i == m, x > 0) %>% ## i = metabolite, j = reaction, x = flux
        dplyr::pull(j)

      r2m_matrix[reactions, m] <- r2m_matrix[reactions,m] + 1
    }

  }else{
    ## Metabolite X maps to reactions where X can be either the reactant or the product
    for(m in 1:ncol(r2m_matrix)){
      reactions <- S %>%
        dplyr::filter(i == m) %>% ## i = metabolite, j = reaction, x = flux
        dplyr::pull(j)

      r2m_matrix[reactions, m] <- r2m_matrix[reactions,m] + 1
    }
  }

  ## association list, name of each element = metabolite, each element = reactions associated with this metabolite
  m2r <- as.data.frame(r2m_matrix) %>%
    lapply(., function(x,y){return(row.names(y)[which(x >0)])}, .)

  ## association list, name of each element = reaction, each element = metabolites associated with this reaction
  r2m <- r2m_matrix %>%
    t() %>%
    as.data.frame() %>%
    lapply(., function(x,y){return(row.names(y)[which(x >0)])}, .)

  ## merged metabolites in "r2m"
  r2m_merged <- lapply(r2m, function(x){x <- x %>%
    substr(., 1, nchar(.)-1) %>%
    unique(.)})

  ## merged metabolites in "m2r"
  unique_met <- unique(substr(names(m2r), 1, nchar(names(m2r))-1))
  m2r_merged <- list()
  for(i in 1:length(unique_met)){
    met <- unique_met[i]
    m2r_merged[[met]] <- m2r[grepl(met, names(m2r))] %>%
      unlist(.) %>%
      unique(.)
  }

  ## map reactions to genes, row = reactions, columns = genes
  data$rxnGeneMat@Dimnames <-  list(reaction_df$name, gene_df$name)
  r2g_matrix <- as(data$rxnGeneMat,"matrix")

  ## association list, name of each element = reaction, each element = genes associated with this reaction
  g2r <- as.data.frame(r2g_matrix) %>%
    lapply(., function(x,y){return(row.names(y)[which(x >0)])}, .)

  ## association list, name of each element = gene ensemble ID, each element = reactions associated with this gene
  r2g <- r2g_matrix %>%
    t(.) %>%
    as.data.frame(.) %>%
    lapply(., function(x,y){return(row.names(y)[which(x >0)])}, .)

  ## association list name of each element = unmerged metabolite, each element = genes associated with this metabolite
  m2g <- list()
  for(i in 1:length(m2r)){
    m2g[[names(m2r)[i]]] <- r2g[m2r[[i]]] %>%
      unlist(.) %>%
      unique(.)
  }

  ## association list name of each element = merged metabolite, each element = genes associated with this metabolite
  m2g_merged <- list()
  for(i in 1:length(m2r_merged)){
    m2g_merged[[names(m2r_merged)[i]]] <- r2g[m2r_merged[[i]]] %>%
      unlist(.) %>%
      unique(.)
  }

  ## association list name of each element = gene, each element = unmerged metabolites associated with this genes
  g2m <- list()
  for(i in 1:length(g2r)){
    g2m[[names(g2r)[i]]] <- r2m[g2r[[i]]] %>%
      unlist(.) %>%
      unique(.)
  }

  ## association list name of each element = gene, each element = merged metabolites associated with this genes
  g2m_merged <- list()
  for(i in 1:length(g2r)){
    g2m_merged[[names(g2r)[i]]] <- r2m_merged[g2r[[i]]] %>%
      unlist(.) %>%
      unique(.)
  }

  return(list(meta_df = meta_df,
              meta_df_merged = NULL,
              reaction_df = reaction_df,
              gene_df = gene_df,
              m2r = m2r,
              m2r_merged =  m2r_merged,
              r2m = r2m,
              r2m_merged = r2m_merged,
              r2g = r2g,
              g2r = g2r,
              m2g = m2g,
              m2g_merged = m2g_merged,
              g2m = g2m,
              g2m_merged = g2m_merged))
}

## intermediate step of mapping metabolites -> reactions, reactions -> genes, and finally, metabolites -> genes
.intermediate  <- function(reactions,g2r){
  mapped_genes <- lapply(g2r[, reactions], function(x){return(which(x!=0, arr.ind = FALSE))}) %>%
    unlist(.) %>%
    unique(.)
  return(mapped_genes)
}
