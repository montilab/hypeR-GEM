#' Title Query and plot the sub-graph
#'
#' @param ig the igraph object derived from a given GEM
#' @param metabolite_query a character vector contains metabolites which are consistent with V(ig)$symbol or V(ig)$name, for metabolite, the "symbol" should be consistent with "fullname" column in hypeR_GEM_obj$mapped_metabolite_signatures
#' @param gene_query a character vector contains enzyme-coding genes which are consistent with V(ig)$symbol or V(ig)$name
#' @param metabolite_by type of symbols used by metabolite_query
#' @param gene_by type of symbols used by gene_query
#' @param focus focus of the subgraph query
#' @param edge_label An edge attribute in E(ig) that can be used to represent edge labels
#' @param main the title of the network
#' @param metabolite_query_color the node color of "metabolite_query"
#' @param gene_query_color the node color of "gene_query"
#' @param associated_metabolite_color the node color of metabolites other than "metabolite_query"
#' @param associated_gene_color the node color of genes other than "gene_query"
#' @param metabolite_shape the node shape of metabolites
#' @param gene_shape the node shape of genes
#' @param hover logical variable
#' @param stabilization logical variable
#' @param igraph_layout the layout algorithm available in "igraph", see https://igraph.org/r/html/1.3.0/layout_.html
#' @param layout_seed random seed of the network layout

#' @import methods utils
#' @importFrom magrittr %>% set_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter rename pull mutate case_when
#' @importFrom igraph as_data_frame neighborhood induced_subgraph vertex_attr edge_attr
#' @importFrom visNetwork visNetwork visNodes visEdges visGroups visOptions visLegend visInteraction visLayout visPhysics visIgraphLayout
#' @importFrom rlang .data
#'
#' @return a list of: ig = node-induced sub-graph, p = a visNetwork object
#' @export
visNet <- function(ig,
                   metabolite_query = NULL,
                   gene_query = NULL,
                   metabolite_by = c("symbol","refmet_name", "HMDB_ID","name"),
                   gene_by = c("symbol","ensemble_ID", "name"),
                   focus = c("all","metabolite","pathway"),
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
  metabolite_by <- match.arg(metabolite_by)
  gene_by <- match.arg(gene_by)
  focus <- match.arg(focus)

  ## check
  if(!(is(ig, "igraph"))) stop(" 'ig' must be a igraph object!\n")
  if(!(is.character(metabolite_query))) stop(" 'metabolite_query' must be a character vector!\n")
  if(!(is.character(gene_query))) stop(" 'gene_query' must be a character vector!\n")
  if(!is.null(edge_label)){ if(is.null(igraph::edge_attr(ig, edge_label))) stop(" 'edge_label' must be an edge attribute in ig!\n") }

  ## query  nodes from "ig"
  df <- igraph::as_data_frame(ig, what='both')

  ## query nodes from "ig"，duplicated "symbol" implies multiple compartments
  query_df <- df$vertices %>%
    magrittr::set_rownames(1:nrow(.data)) %>%
    tibble::rownames_to_column(var='vid') %>%
    dplyr::mutate(vid = as.numeric(.data$vid)) %>%
    dplyr::filter(!is.na(!!as.name(.data$metabolite_by)) & !is.na(!!as.name(.data$gene_by))) %>%
    dplyr::filter(!!as.name(.data$metabolite_by) %in% metabolite_query | !!as.name(.data$gene_by) %in% gene_query)


  ## obtain node Ids of gene_query，empty if gene_query = NULL
  gene_query_vids <- query_df %>%
    dplyr::filter(.data$node_type == "gene") %>%
    dplyr::pull(.data$vid)

  ## obtain node Ids of metabolite_query，empty if metabolite_query = NULL
  metabolite_query_vids <- query_df %>%
    dplyr::filter(.data$node_type == "metabolite") %>%
    dplyr::pull(.data$vid)

  ## gene_query and their neighbor, empty if gene_query = NULL
  gene_neighbor_vids <- lapply(igraph::neighborhood(ig, order=1, nodes= gene_query_vids, mode='all'), as.numeric) %>%
    unlist() %>%
    unique()

  ## metabolite_query and their neighbor，empty if metabolite_query = NULL
  metabolite_neighbor_vids <- lapply(igraph::neighborhood(ig, order=1, nodes= metabolite_query_vids, mode='all'), as.numeric) %>%
    unlist() %>%
    unique()

  ## specified visualization focus
  if(focus == "all"){
    all_vids <- unique(c(gene_neighbor_vids, metabolite_neighbor_vids))

  }
  if(focus == "metabolite"){
    all_vids <-  unique(c(metabolite_neighbor_vids, gene_query_vids))
  }
  if(focus == "pathway"){
    all_vids <- unique(c(gene_neighbor_vids, metabolite_query_vids))


  }

  ## Obtain the node-induced sub-graph
  ig_sub <- igraph::induced_subgraph(ig, vids= all_vids)

  if(!is.null(vertex_attr(ig_sub, 'id'))){
    nodes <- igraph::as_data_frame(ig_sub, what = 'vertices')
    nodes[nodes$symbol %in% metabolite_query, 'node_type'] <- 'metabolite_query'
    nodes[nodes$symbol %in% gene_query, 'node_type'] <- 'gene_query'

    nodes <- nodes %>%
      dplyr::mutate(compartment = case_when(
        .data$node_type %in% c("metabolite","metabolite_query") ~ stringr::str_sub(.data$name, -1),
        TRUE ~ "")) %>%
      dplyr::mutate(symbol = case_when(
        .data$node_type %in% c("metabolite","metabolite_query") ~ paste(.data$symbol, .data$compartment, sep = "_"),
        TRUE ~ symbol
      )) %>%
      dplyr::rename(group = .data$node_type,
                    label = .data$symbol,
                    old_id = .data$id) %>%
      dplyr::filter(!is.na(.data$label)) %>%
      tibble::rownames_to_column(var='id')

  }else{
    nodes <- igraph::as_data_frame(ig_sub, what = 'vertices')
    nodes[nodes$symbol %in% metabolite_query, 'node_type'] <- 'metabolite_query'
    nodes[nodes$symbol %in% gene_query, 'node_type'] <- 'gene_query'

    nodes <- nodes %>%
      dplyr::mutate(compartment = case_when(
        .data$node_type %in% c("metabolite","metabolite_query") ~ stringr::str_sub(.data$name, -1),
        TRUE ~ "")) %>%
      dplyr::mutate(symbol = case_when(
        .data$node_type %in% c("metabolite","metabolite_query") ~ paste(.data$symbol, .data$compartment, sep = "_"),
        TRUE ~ .data$symbol
      )) %>%
      dplyr::rename(group = .data$node_type,
                    label = .data$symbol) %>%
      dplyr::filter(!is.na(.data$label)) %>%
      tibble::rownames_to_column(var='id')

  }

  ## create "edge" dataframe
  if(!is.null(edge_label)){
    edges <- igraph::as_data_frame(ig_sub, what='edges') %>%
      dplyr::rename(label = !!as.name(.data$edge_label))
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
              nodes = nodes,
              edges = edges,
              p=p))

}
