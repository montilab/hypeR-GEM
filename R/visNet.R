#' Title Query and plot the sub-graph
#'
#' @param ig the igraph object derived from a given GEM
#' @param metabolite_query a character vector contains metabolites which are consistent with V(ig)$symbol or V(ig)$name, for metabolite, the "symbol" should be consistent with "fullname" column in hypeR_GEM_obj$mapped_metabolite_signatures
#' @param gene_query a character vector contains enzyme-coding genes which are consistent with V(ig)$symbol or V(ig)$name
#' @param by if the queries are given by V(ig)$symbol, then by = "symbol"; if the queries are given by V(ig)$name, then by = "name"
#' @param edge_label An edge attribute in E(ig) that can be used to represent edge labels
#' @param main the title of the network
#' @param metabolite_query_color the node color of "metabolite_query"
#' @param gene_query_color the node color of "gene_query"
#' @param associated_metabolite_color the node color of metabolites other than "metabolite_query"
#' @param associated_gene_color the node color of genes other than "gene_query"
#' @param metabolite_shape the node shape of metabolites
#' @param gene_shape the node shape of genes
#' @param igraph_layout the layout algorithm available in "igraph", see https://igraph.org/r/html/1.3.0/layout_.html
#' @param layout_seed random seed of the network layout

#'@importFrom magrittr %>%
#'@importFrom tibble rownames_to_column
#'@importFrom dplyr filter rename pull
#'@importFrom igraph as_data_frame neighborhood induced_subgraph
#'@importFrom visNetwork visNetwork visNodes visEdges visGroups visOptions visLegend visInteraction visLayout visPhysics visIgraphLayout

#' @return a list of: ig = node-induced sub-graph, p = a visNetwork object
#' @export
visNet <- function(ig,
                   metabolite_query = NULL,
                   gene_query = NULL,
                   by = c("symbol","name"),
                   edge_label = NULL,
                   main = "A simple reaction-based network",
                   metabolite_query_color = "#FF9999",
                   gene_query_color = "#56B4E9",
                   associated_metabolite_color = "#808080",
                   associated_gene_color = "#808080",
                   metabolite_shape = 'dot',
                   gene_shape = 'diamond',
                   hover = TRUE,
                   stabilization = FALSE,
                   igraph_layout = "layout_with_fr",
                   layout_seed = 42){

  # Default arguments
  by <- match.arg(by)

  ## query  nodes from "ig"
  df <- igraph::as_data_frame(ig,what='both')

  ## query nodes from "ig"，duplicated "symbol" implies multiple compartments
  query_df <- df$vertices %>%
    tibble::rownames_to_column(var='vid') %>%
    dplyr::filter(!is.na(!!as.name(by))) %>%
    dplyr::filter(!!as.name(by) %in% metabolite_query | !!as.name(by) %in% gene_query)

  ## obtain neighbors of gene_query
  gene_query_vids <- query_df %>%
    dplyr::filter(node_type == "gene") %>%
    dplyr::pull(vid)

  metabolite_query_vids <- query_df %>%
    dplyr::filter(node_type == "metabolite") %>%
    dplyr::pull(vid)

  ## metabolites belong to"metabolite_query" and other metabolites
  gene_neighbor_vids <- lapply(igraph::neighborhood(ig, order=1, nodes=gene_query_vids, mode='all'), as.numeric) %>%
    unlist(.) %>%
    unique(.)

  ## metabolite_query and their neighbor
  metabolite_neighbor_vids <- lapply(igraph::neighborhood(ig, order=1, nodes=metabolite_query_vids, mode='all'), as.numeric) %>%
    unlist(.) %>%
    unique(.)

  all_neighbor_vids <- unique(c(gene_neighbor_vids, metabolite_neighbor_vids))

  ## Obtain the node-induced sub-graph
  ig_sub <- igraph::induced_subgraph(ig, vids= all_neighbor_vids)

  if(!is.null(vertex_attr(ig_sub, 'id'))){
    nodes <- igraph::as_data_frame(ig_sub, what='vertices') %>%
      dplyr::rename(group = node_type,
                    label = symbol,
                    old_id = id) %>%
      dplyr::filter(!is.na(label)) %>%
      tibble::rownames_to_column(var='id')

    nodes[nodes$label %in% metabolite_query, 'group'] <- 'metabolite_query'
    nodes[nodes$label %in% gene_query, 'group'] <- 'gene_query'
  }else{
    nodes <- igraph::as_data_frame(ig_sub, what='vertices') %>%
      dplyr::rename(group = node_type,
                    label = symbol) %>%
      dplyr::filter(!is.na(label)) %>%
      tibble::rownames_to_column(var='id')
    nodes[nodes$label %in% metabolite_query, 'group'] <- 'metabolite_query'
    nodes[nodes$label %in% gene_query, 'group'] <- 'gene_query'
  }

  ## create "edge" dataframe
  if(!is.null(edge_label)){
    edges <- igraph::as_data_frame(ig_sub, what='edges') %>%
      dplyr::rename(label = !!as.name(edge_label))
  }else{
    edges <- igraph::as_data_frame(ig_sub, what='edges')
  }

  ## Add values to nodes
  p <- visNetwork::visNetwork(nodes, edges, main=main) %>%
    visNetwork::visNodes(id='id',label='label', size=25) %>%
    visNetwork::visEdges(arrows ="to") %>%
    visNetwork::visGroups(groupname = "metabolite_query", color = metabolite_query_color, shape = metabolite_shape) %>%
    visNetwork::visGroups(groupname = "gene_query", color = gene_query_color, shape = gene_shape) %>%
    visNetwork::visGroups(groupname = "metabolite", color = associated_metabolite_color, shape = metabolite_shape) %>%
    visNetwork::visGroups(groupname = "gene", color = associated_gene_color, shape = gene_shape) %>%
    visNetwork::visOptions(highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visLegend() %>%
    visNetwork::visInteraction(hover = hover) %>%
    visNetwork::visPhysics(stabilization = stabilization) %>%
    visNetwork::visIgraphLayout(layout = igraph_layout) %>%
    visNetwork::visLayout(randomSeed = layout_seed)

  return(list(ig=ig_sub,
              p=p))

}
