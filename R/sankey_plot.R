#' Create a 3-layer Sankey diagram (Metabolites → MSETs → GSETs)
#'
#' Builds a Sankey diagram with three layers based on a hypeR-GEM object:
#' metabolites → metabolite sets (MSETs) → gene sets (GSETs).
#'
#' @param hypeR_GEM_obj A list with components \code{$mapped_metabolite_signatures}
#'   (a list of data frames) and \code{$gene_tables} (a list of data frames).
#' @param msets_list A list of metabolite-set collections. Each element should be
#'   a named list where names are set names and values are character vectors of
#'   metabolite identifiers matching \code{key}.
#' @param gsets_list A list of gene-set collections. Each element should be a
#'   named list where names are set names and values are character vectors of genes.
#' @param key Column name in your metabolite tables used as the metabolite
#'   identifier (e.g., \code{"refmet_name"}, \code{"HMDB"}). Default: \code{"refmet_name"}.
#' @param font_size Numeric font size for node labels. Default: 12.
#' @param node_width Numeric node width. Default: 30.
#'
#' @return A \code{networkD3::sankeyNetwork} htmlwidget.
#' @export
#'
#' @importFrom networkD3 sankeyNetwork
#' @importFrom dplyr bind_rows select mutate filter group_by summarise rename right_join
#'   distinct pull all_of any_of
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_trim str_c
#' @importFrom rlang sym .data
#' @importFrom stats na.omit
sankey_plot <- function(hypeR_GEM_obj,
                        msets_list,
                        gsets_list,
                        key = "refmet_name",
                        font_size = 12,
                        node_width = 30) {

  # ---------- basic checks ----------
  if (!is.list(hypeR_GEM_obj) ||
      is.null(hypeR_GEM_obj$mapped_metabolite_signatures) ||
      is.null(hypeR_GEM_obj$gene_tables)) {
    stop("`hypeR_GEM_obj` must contain $mapped_metabolite_signatures and $gene_tables.")
  }
  if (!is.list(msets_list) || !length(msets_list)) {
    stop("`msets_list` must be a non-empty list of MSET collections.")
  }
  if (!is.list(gsets_list) || !length(gsets_list)) {
    stop("`gsets_list` must be a non-empty list of GSET collections.")
  }
  if (!is.character(key) || length(key) != 1L) {
    stop("`key` must be a single character column name.")
  }

  # remove empty inner elements
  msets_list <- msets_list[lengths(msets_list) > 0L]
  gsets_list <- gsets_list[lengths(gsets_list) > 0L]

  # ---------- build edges and tables ----------
  met2MSET_edge_list <- .met2MSET_links(hypeR_GEM_obj, msets_list, key = key)

  met_table   <- .create_met_table(hypeR_GEM_obj, msets_list, key = key)
  gene_table  <- .create_gene_table(hypeR_GEM_obj, gsets_list)
  msets_table <- .create_msets_table(met_table, key = key)
  gsets_table <- .create_gsets_table(gene_table)

  msets2gsets_edge_list <- .msets2gsets_links(
    met_table, msets_table, gene_table, gsets_table, key = key
  )

  links <- rbind(met2MSET_edge_list, msets2gsets_edge_list)
  if (!nrow(links)) stop("No edges found to plot (all overlaps were zero).")

  # ---------- node index ----------
  nodes <- data.frame(name = union(unique(links$source), unique(links$target)),
                      stringsAsFactors = FALSE)
  links$source <- match(links$source, nodes$name) - 1L
  links$target <- match(links$target, nodes$name) - 1L

  if (any(is.na(links$source)) || any(is.na(links$target)) ||
      any(links$source < 0L) || any(links$target < 0L)) {
    stop("Internal error: node indexing failed.")
  }

  networkD3::sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "source", Target = "target", Value = "weight",
    NodeID = "name", fontSize = font_size, nodeWidth = node_width,
    sinksRight = FALSE
  )
}

# ----- internals -----------------------------------------------------------

#' @keywords internal
.met2MSET_links <- function(hypeR_GEM_obj, msets_list, key = "refmet_name") {
  met_df <- dplyr::bind_rows(hypeR_GEM_obj$mapped_metabolite_signatures,
                             .id = "signature_type")

  if (!all(c("signature_type", key) %in% names(met_df))) {
    stop("`mapped_metabolite_signatures` must include columns 'signature_type' and '", key, "'.")
  }

  # split metabolites by signature_type
  met_groups <- split(met_df[[key]], met_df[["signature_type"]])
  met_groups <- lapply(met_groups, function(v) unique(stats::na.omit(as.character(v))))

  # flatten one level: list(named vectors)
  msets <- unlist(msets_list, recursive = FALSE, use.names = TRUE)
  msets <- lapply(msets, function(v) unique(stats::na.omit(as.character(v))))
  if (is.null(names(msets))) names(msets) <- paste0("MSET_", seq_along(msets))

  edge_list <- base::expand.grid(
    source = base::names(met_groups),
    target = base::names(msets),
    stringsAsFactors = FALSE
  )
  edge_list$weight <- mapply(function(i, j) {
    base::length(base::intersect(met_groups[[i]], msets[[j]]))
  }, edge_list$source, edge_list$target)

  edge_list[edge_list$weight > 0L, , drop = FALSE]
}

#' @keywords internal
.find_association <- function(id, sets) {
  sets <- lapply(sets, function(x) unlist(x, use.names = FALSE))
  is_member <- vapply(sets, function(x) id %in% x, logical(1))
  names(sets)[is_member]
}

#' @keywords internal
.find_associated_sets <- function(id, sets_list) {
  out <- unlist(lapply(sets_list, .find_association, id = id), use.names = TRUE)
  paste(out, collapse = ";")
}

#' @keywords internal
.create_met_table <- function(hypeR_GEM_obj, msets_list, key = "refmet_name", sep = ";") {
  k <- rlang::sym(key)

  # metabolites + associated MSETs (by signature type)
  met_df <- dplyr::bind_rows(hypeR_GEM_obj$mapped_metabolite_signatures, .id = "type") |>
    dplyr::select(!!k, dplyr::all_of("type")) |>
    dplyr::mutate(
      associated_msets = vapply(
        X = .data[[key]],
        FUN = .find_associated_sets,
        FUN.VALUE = character(1),
        sets_list = msets_list
      )
    )

  gene_tab <- dplyr::bind_rows(hypeR_GEM_obj$gene_tables, .id = "type")
  if (!all(c("symbol", "associated_metabolites") %in% names(gene_tab))) {
    stop("`gene_tables` must contain columns 'symbol' and 'associated_metabolites'.")
  }

  met_table <- gene_tab |>
    dplyr::select(dplyr::all_of(c("symbol", "associated_metabolites"))) |>
    tidyr::separate_rows(dplyr::all_of("associated_metabolites"), sep = sep) |>
    dplyr::mutate(associated_metabolites = stringr::str_trim(.data$associated_metabolites)) |>
    dplyr::filter(.data$associated_metabolites != "", !is.na(.data$associated_metabolites)) |>
    dplyr::group_by(metabolite = .data$associated_metabolites) |>
    dplyr::summarise(
      associated_genes = paste(unique(.data$symbol), collapse = sep),
      .groups = "drop"
    ) |>
    dplyr::rename(!!k := rlang::sym("metabolite")) |>
    dplyr::right_join(met_df, by = key)

  met_table
}

#' @keywords internal
.create_msets_table <- function(met_table, key = "refmet_name") {
  k <- rlang::sym(key)

  met_table |>
    dplyr::filter(.data$associated_msets != "") |>
    dplyr::select(!!k, dplyr::all_of("associated_msets")) |>
    tidyr::separate_rows(dplyr::all_of("associated_msets"), sep = ";") |>
    dplyr::filter(nzchar(.data$associated_msets)) |>
    dplyr::group_by(.data$associated_msets) |>
    dplyr::summarise(
      associated_metabolite = stringr::str_c(.data[[key]], collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::rename(msets = rlang::sym("associated_msets"))
}

#' @keywords internal
.create_gene_table <- function(hypeR_GEM_obj, gsets_list) {
  gene_df <- dplyr::bind_rows(hypeR_GEM_obj$gene_tables, .id = "gene_type") |>
    dplyr::select(dplyr::any_of(c("name", "symbol")))

  if (!("symbol" %in% names(gene_df))) {
    stop("`gene_tables` must contain a 'symbol' column.")
  }

  gene_df$associated_gsets <- vapply(
    X = gene_df$symbol,
    FUN = .find_associated_sets,
    FUN.VALUE = character(1),
    sets_list = gsets_list
  )
  gene_df
}

#' @keywords internal
.create_gsets_table <- function(gene_table) {
  if (!all(c("symbol", "associated_gsets") %in% names(gene_table))) {
    stop("`gene_table` must contain 'symbol' and 'associated_gsets'.")
  }

  gene_table |>
    dplyr::select(dplyr::all_of(c("symbol", "associated_gsets"))) |>
    dplyr::filter(.data$associated_gsets != "") |>
    tidyr::separate_rows(dplyr::all_of("associated_gsets"), sep = ";") |>
    dplyr::group_by(.data$associated_gsets) |>
    dplyr::summarise(
      associated_gene = stringr::str_c(.data$symbol, collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::rename(gsets = rlang::sym("associated_gsets"))
}

#' @keywords internal
.msets2gsets_links <- function(met_table, msets_table, gene_table, gsets_table, key = "refmet_name") {
  edge_list <- data.frame(source = character(), target = character(),
                          weight = integer(), stringsAsFactors = FALSE)

  if (!all(c("msets", "associated_metabolite") %in% names(msets_table))) return(edge_list)
  if (!all(c(key, "associated_genes") %in% names(met_table))) return(edge_list)
  if (!all(c("symbol", "associated_gsets") %in% names(gene_table))) return(edge_list)
  if (!all(c("gsets", "associated_gene") %in% names(gsets_table))) return(edge_list)

  # Precompute: gset -> gene vector
  gset_genes <- lapply(seq_len(nrow(gsets_table)), function(i) {
    gs <- unique(unlist(strsplit(stats::na.omit(gsets_table$associated_gene[i]),
                                 ";", fixed = TRUE), use.names = FALSE))
    gs[nzchar(gs)]
  })
  names(gset_genes) <- as.character(gsets_table$gsets)

  for (i in seq_len(nrow(msets_table))) {
    mset_i <- as.character(msets_table$msets[i])
    assoc_met_str <- as.character(msets_table$associated_metabolite[i])
    if (!nzchar(mset_i) || !nzchar(assoc_met_str)) next

    mets <- unlist(strsplit(assoc_met_str, ";", fixed = TRUE), use.names = FALSE)
    mets <- mets[nzchar(mets)]
    if (!length(mets)) next

    # genes from these metabolites
    mt_sub <- met_table[met_table[[key]] %in% mets, , drop = FALSE]
    ag <- unlist(strsplit(stats::na.omit(mt_sub$associated_genes), ";", fixed = TRUE),
                 use.names = FALSE)
    m_genes <- unique(ag[nzchar(ag)])
    if (!length(m_genes)) next

    # which gsets these genes belong to (by gene_table)
    gt_sub <- gene_table[gene_table$symbol %in% m_genes, , drop = FALSE]
    gs_flat <- unlist(strsplit(stats::na.omit(gt_sub$associated_gsets), ";", fixed = TRUE),
                      use.names = FALSE)
    gsets <- unique(gs_flat[nzchar(gs_flat)])
    if (!length(gsets)) next

    w <- vapply(gsets, function(g) {
      gs_genes <- gset_genes[[g]]
      if (is.null(gs_genes)) 0L else length(intersect(m_genes, gs_genes))
    }, integer(1))

    keep <- w > 0L
    if (any(keep)) {
      edge_list <- rbind(
        edge_list,
        data.frame(source = rep(mset_i, sum(keep)),
                   target = gsets[keep],
                   weight = as.integer(w[keep]),
                   stringsAsFactors = FALSE)
      )
    }
  }

  dplyr::distinct(edge_list)
}

# Silence R CMD check for data-masked column names used in verbs
utils::globalVariables(c(
  "type", "symbol", "associated_metabolites", "metabolite",
  "associated_msets", "associated_gsets", "name"
))
