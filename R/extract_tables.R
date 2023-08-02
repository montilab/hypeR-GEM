#' Title extract tables from a given genome scale metabolic model(GEM) of ".mat" format
#'
#' @param con a connection or the path of the GEM
#' @return a list of tables

#'@importFrom R.matlab readMat
#'@importFrom magrittr %>%
#'@importFrom Matrix Matrix summary
#'@importFrom dplyr arrange filter select pull

#' @examples
#' con = file.path("data/Recon3D.mat")
#' tables = extract_tables(con)

#' @export
extract_tables <- function(con){

  ## load GEM and extract data
  GEM <- R.matlab::readMat(con=file.path(con))
  data <- GEM[[names(GEM)]][,,1]

  ## create basic data frames of metabolites, reactions, and genes
  compos <- sapply(data$comps[,1],unlist)
  meta_df <- data.frame(name = sapply(data$mets[,1],unlist),
                        fullname = sapply(data$metNames[,1],unlist),
                        compartment = compos[sapply(data$metComps[,1],unlist)])

  reaction_df <- data.frame(name = sapply(data$rxns[,1],unlist))

  gene_df <- data.frame(name = unlist(sapply(data$genes,unlist)))

  ## map metabolites to reactions, row = reactions, column = metabolites
  r2m <- matrix(0,nrow=nrow(reaction_df),
                ncol=nrow(meta_df),
                dimnames = list(reaction_df$name,meta_df$name))

  ### pull out the S matrix, Matrix::summary() converts S@i, S@j from 0-based to 1-based
  S <- as(as(as(data$S, "dMatrix"), "generalMatrix"), "TsparseMatrix") %>%
    Matrix::summary(.) %>%
    as.data.frame(.) %>%
    dplyr::arrange(i)

  for(m in 1:ncol(r2m)){
    reactions <- S %>%
      dplyr::filter(i == m) %>%
      dplyr::pull(j)

    r2m[reactions, m] <- r2m[reactions,m] + 1
  }
  r2m <- as.data.frame(r2m)
  #r2m <- Matrix::Matrix(r2m, sparse = TRUE)

  ## map reactions to genes, row = genes, columns = reactions
  data$rxnGeneMat@Dimnames <-  list(reaction_df$name, gene_df$name)
  g2r <- as(data$rxnGeneMat,"matrix") %>%
    t(.) %>%
    as.data.frame(.)

  ## map genes to metabolites, row = metabolites, column = genes
  m2g <- matrix(0,nrow=nrow(meta_df),
                ncol=nrow(gene_df),
                dimnames = list(meta_df$name,gene_df$name))

  m2g_list <- lapply(r2m, function(x){return(which(x!=0, arr.ind = FALSE))}) %>%
    lapply(., .intermediate, g2r)

  for(m in 1:length(m2g_list)){
    genes <- m2g_list[[m]]
    m2g[m, genes] <- m2g[m, genes] + 1
  }
  #m2g <- Matrix::Matrix(m2g, sparse = TRUE)

  return(list(meta_df = meta_df,
              reaction_df = reaction_df,
              gene_df = gene_df,
              r2m = Matrix::Matrix(as(r2m,"matrix"), sparse = TRUE),
              g2r = Matrix::Matrix(as(g2r,"matrix"), sparse = TRUE),
              m2g = Matrix::Matrix(m2g, sparse = TRUE)))

}

## intermediate step of mapping metabolites -> reactions, reactions -> genes, and finally, metabolites -> genes
.intermediate  <- function(reactions,g2r){
  mapped_genes <- lapply(g2r[, reactions], function(x){return(which(x!=0, arr.ind = FALSE))}) %>%
    unlist(.) %>%
    unique(.)
  return(mapped_genes)
}
