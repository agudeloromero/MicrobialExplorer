# =============================================================================
# MICROBIOME NETWORK ANALYSIS
# =============================================================================
# Description : Comprehensive co-occurrence network analysis pipeline covering
#               network construction, topology, hub identification, community
#               detection, differential networks, and random network comparison
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Network plots, topology tables, hub taxa, community modules,
#               differential network comparisons
# Author      : Patricia
# Dependencies: phyloseq, igraph, ggraph, tidygraph, ggplot2, dplyr, tidyr,
#               patchwork, RColorBrewer, vegan, scales, tibble, stringr,
#               SPIEC-EASI (optional), SpiecEasi (optional)
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(vegan)
  library(scales)
  library(tibble)
  library(stringr)
  library(forcats)
})

pkg_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Global theme
theme_microbiome <- function() {
  theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#2c3e50", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold", size = 10),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      legend.title      = element_text(size = 10, face = "bold"),
      legend.text       = element_text(size = 9),
      plot.title        = element_text(size = 13, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = 10, colour = "grey40", hjust = 0),
      plot.caption      = element_text(size = 8, colour = "grey60", hjust = 1)
    )
}

theme_network <- function() {
  theme_void() +
    theme(
      legend.title      = element_text(size = 10, face = "bold"),
      legend.text       = element_text(size = 9),
      plot.title        = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = 10, colour = "grey40", hjust = 0.5),
      plot.caption      = element_text(size = 8, colour = "grey60", hjust = 1),
      plot.background   = element_rect(fill = "white", colour = NA)
    )
}


# =============================================================================
# SECTION 1 — DATA PREPARATION FOR NETWORK ANALYSIS
# =============================================================================

#' Prepare a phyloseq object for network analysis.
#'
#' Filters to prevalent taxa, applies CLR transformation appropriate
#' for compositional network inference, and optionally subsets samples.
#'
#' @param ps              A phyloseq object (raw counts).
#' @param rank            Taxonomic rank to agglomerate to. Default = "Genus".
#' @param min_prevalence  Minimum fraction of samples a taxon must appear in.
#' @param min_abundance   Minimum mean relative abundance (0–1). Default = 0.001.
#' @param max_taxa        Maximum taxa to retain (by prevalence). Default = 150.
#' @param group_var       Optional group variable for subsetting.
#' @param group_subset    Optional group level to subset to.
#' @return A filtered phyloseq object with CLR-transformed matrix attached.

prepare_network_data <- function(ps,
                                  rank           = "Genus",
                                  min_prevalence = 0.20,
                                  min_abundance  = 0.001,
                                  max_taxa       = 150,
                                  group_var      = NULL,
                                  group_subset   = NULL) {

  cat("=== Preparing data for network analysis ===\n")

  # Optional group subsetting
  if (!is.null(group_var) && !is.null(group_subset)) {
    ps <- subset_samples(ps,
                          as.character(sample_data(ps)[[group_var]]) == group_subset)
    cat("  Subset to", group_subset, ":", nsamples(ps), "samples\n")
  }

  # Agglomerate to rank
  ps_agg <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  dup_idx <- duplicated(new_names)
  new_names <- make.unique(new_names, sep = "_dup")
  taxa_names(ps_agg) <- new_names

  # Relative abundance and filter
  ps_rel  <- transform_sample_counts(ps_agg, function(x) x / sum(x))
  otu_mat <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) otu_mat <- t(otu_mat)

  prevalence  <- rowSums(otu_mat > 0) / ncol(otu_mat)
  mean_abund  <- rowMeans(otu_mat)
  keep        <- prevalence >= min_prevalence & mean_abund >= min_abundance

  # Further trim to max_taxa by prevalence
  if (sum(keep) > max_taxa) {
    keep_idx <- which(keep)
    top_idx  <- keep_idx[order(prevalence[keep_idx], decreasing = TRUE)][seq_len(max_taxa)]
    keep     <- seq_len(nrow(otu_mat)) %in% top_idx
  }

  ps_filt <- prune_taxa(keep, ps_agg)
  cat("  Taxa after filtering:", ntaxa(ps_filt), "\n")
  cat("  Samples:", nsamples(ps_filt), "\n\n")

  return(ps_filt)
}


# =============================================================================
# SECTION 2 — NETWORK CONSTRUCTION
# =============================================================================

#' Build a co-occurrence network using Spearman correlation with
#' multiple testing correction and optional SPIEC-EASI support.
#'
#' @param ps              A phyloseq object.
#' @param method          "spearman", "pearson", or "spiec-easi". Default = "spearman".
#' @param cor_threshold   Minimum absolute correlation to include an edge. Default = 0.6.
#' @param p_threshold     Maximum adjusted p-value for an edge. Default = 0.05.
#' @param p_adjust_method Multiple testing correction. Default = "BH".
#' @param seed            Random seed. Default = 42.
#' @return A list: igraph object, correlation matrix, and adjacency matrix.

build_network <- function(ps,
                           method          = "spearman",
                           cor_threshold   = 0.3,
                           p_threshold     = 0.05,
                           p_adjust_method = "BH",
                           seed            = 42) {

  cat("=== Building co-occurrence network ===\n")
  cat("  Method:", method, "\n")

  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

  # CLR transformation for compositional data
  otu_clr <- t(apply(otu_mat + 0.5, 2, function(x) {
    log(x) - mean(log(x))
  }))   # samples as rows

  n_taxa <- nrow(otu_mat)

  if (method == "spiec-easi") {
    # --- SPIEC-EASI (optimal for compositional data) -----------------------
    if (!pkg_available("SpiecEasi")) {
      cat("  SpiecEasi not installed - falling back to Spearman\n")
      cat("  Install with: devtools::install_github('zdk123/SpiecEasi')\n")
      method <- "spearman"
    } else {
      library(SpiecEasi)
      cat("  Running SPIEC-EASI (MB method)...\n")
      set.seed(seed)
      se_res   <- spiec.easi(ps, method = "mb",
                              lambda.min.ratio = 1e-2,
                              nlambda = 20,
                              pulsar.params = list(rep.num = 50,
                                                   ncores  = 1))
      adj_mat  <- as.matrix(getRefit(se_res))
      cor_mat  <- as.matrix(cov2cor(as.matrix(getOptCov(se_res))))
      rownames(adj_mat) <- colnames(adj_mat) <- taxa_names(ps)
      rownames(cor_mat) <- colnames(cor_mat) <- taxa_names(ps)
      diag(adj_mat)     <- 0
    }
  }

  if (method %in% c("spearman", "pearson")) {
    # --- Correlation-based network ----------------------------------------
    cat("  Computing pairwise", method, "correlations...\n")

    # Correlation matrix
    cor_mat <- cor(otu_clr, method = method)

    # Significance testing for each pair
    n_samples <- nrow(otu_clr)
    p_mat     <- matrix(1, nrow = n_taxa, ncol = n_taxa)
    rownames(p_mat) <- colnames(p_mat) <- taxa_names(ps)

    for (i in seq_len(n_taxa - 1)) {
      for (j in seq(i + 1, n_taxa)) {
        ct <- cor.test(otu_clr[, i], otu_clr[, j],
                       method = method, exact = FALSE)
        p_mat[i, j] <- ct$p.value
        p_mat[j, i] <- ct$p.value
      }
    }

    # BH correction (only upper triangle)
    upper_idx       <- upper.tri(p_mat)
    p_adj           <- p_mat
    p_adj[upper_idx] <- p.adjust(p_mat[upper_idx], method = p_adjust_method)
    p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]

    # Adjacency: keep edges where |r| >= threshold AND p.adj < threshold
    adj_mat <- matrix(0, nrow = n_taxa, ncol = n_taxa)
    rownames(adj_mat) <- colnames(adj_mat) <- taxa_names(ps)
    sig_mask <- abs(cor_mat) >= cor_threshold & p_adj < p_threshold
    diag(sig_mask) <- FALSE
    adj_mat[sig_mask] <- cor_mat[sig_mask]
  }

  # --- Build igraph ---------------------------------------------------------
  cat("  Building igraph object...\n")

  # Separate positive and negative edges
  pos_adj <- adj_mat
  neg_adj <- adj_mat
  pos_adj[pos_adj < 0] <- 0
  neg_adj[neg_adj > 0] <- 0

  g <- graph_from_adjacency_matrix(
    abs(adj_mat),
    mode     = "undirected",
    weighted = TRUE,
    diag     = FALSE
  )

  # Add edge sign attribute
  edge_list     <- as_edgelist(g)
  edge_signs    <- sapply(seq_len(nrow(edge_list)), function(i) {
    sign(adj_mat[edge_list[i, 1], edge_list[i, 2]])
  })
  E(g)$sign     <- edge_signs
  E(g)$color_e  <- ifelse(edge_signs > 0, "#27ae60", "#e74c3c")

  # Remove isolated nodes
  g <- delete_vertices(g, igraph::degree(g) == 0)

  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  n_pos   <- sum(E(g)$sign > 0)
  n_neg   <- sum(E(g)$sign < 0)

  cat("  Network nodes:", n_nodes, "\n")
  cat("  Network edges:", n_edges, "(", n_pos, "positive,", n_neg, "negative)\n")
  cat("  Network density:", round(igraph::edge_density(g), 4), "\n\n")

  return(list(
    graph      = g,
    cor_mat    = cor_mat,
    adj_mat    = adj_mat,
    method     = method,
    n_nodes    = n_nodes,
    n_edges    = n_edges
  ))
}


# =============================================================================
# SECTION 3 — NETWORK TOPOLOGY
# =============================================================================

#' Calculate comprehensive network topology metrics per node and globally.
#'
#' @param g         An igraph object.
#' @param ps        A phyloseq object (for adding taxonomy to node table).
#' @param rank      Taxonomic rank label. Default = "Genus".
#' @return A list: node metrics data frame and global topology summary.

calculate_network_topology <- function(g, ps = NULL, rank = "Genus") {

  cat("=== Calculating network topology ===\n")

  # --- Global metrics -------------------------------------------------------
  global <- list(
    n_nodes          = vcount(g),
    n_edges          = ecount(g),
    density          = round(igraph::edge_density(g), 4),
    mean_degree      = round(mean(igraph::degree(g)), 3),
    clustering_coef  = round(igraph::transitivity(g, type = "global"), 4),
    avg_path_length  = round(igraph::mean_distance(g, directed = FALSE), 4),
    diameter         = igraph::diameter(g, directed = FALSE),
    modularity       = NA,
    n_components     = igraph::components(g)$no,
    pct_positive     = round(100 * sum(E(g)$sign > 0) / ecount(g), 1)
  )

  # Modularity via Louvain community detection
  set.seed(42)
  comm <- cluster_louvain(g)
  global$modularity <- round(igraph::modularity(comm), 4)
  global$n_modules  <- length(unique(igraph::membership(comm)))

  cat("  Nodes:", global$n_nodes, "| Edges:", global$n_edges, "\n")
  cat("  Density:", global$density, "| Clustering:", global$clustering_coef, "\n")
  cat("  Modularity:", global$modularity, "| Modules:", global$n_modules, "\n\n")

  # --- Node-level metrics ---------------------------------------------------
  node_df <- data.frame(
    taxon           = V(g)$name,
    degree          = igraph::degree(g),
    strength        = strength(g),
    betweenness     = round(igraph::betweenness(g, normalized = TRUE), 6),
    closeness       = round(igraph::closeness(g, normalized = TRUE), 6),
    eigenvector     = round(igraph::eigen_centrality(g)$vector, 6),
    clustering_node = round(igraph::transitivity(g, type = "local"), 4),
    module          = igraph::membership(comm),
    stringsAsFactors = FALSE
  )

  # Hub score: composite of degree, betweenness, and eigenvector
  node_df <- node_df %>%
    mutate(
      degree_z      = scale(degree)[,1],
      between_z     = scale(betweenness)[,1],
      eigen_z       = scale(eigenvector)[,1],
      hub_score     = round((degree_z + between_z + eigen_z) / 3, 4),
      is_hub        = hub_score > quantile(hub_score, 0.9, na.rm = TRUE)
    )

  # Add taxonomy if phyloseq provided
  if (!is.null(ps)) {
    tax_df <- tryCatch({
      raw_mat <- ps@tax_table@.Data
      df <- data.frame(raw_mat, stringsAsFactors = FALSE)
      rownames(df) <- taxa_names(ps)
      colnames(df) <- colnames(ps@tax_table)
      df %>% rownames_to_column("taxon") %>%
        select(taxon, any_of(c("Phylum", "Class", "Family", rank)))
    }, error = function(e) NULL)
    if (!is.null(tax_df)) node_df <- left_join(node_df, tax_df, by = "taxon")
  }


  node_df <- node_df %>% arrange(desc(hub_score))

  n_hubs <- sum(node_df$is_hub, na.rm = TRUE)
  cat("  Hub taxa (top 10% hub score):", n_hubs, "\n")
  cat("  Top 5 hubs:\n")
  print(node_df %>% filter(is_hub) %>%
          select(taxon, degree, betweenness, hub_score) %>%
          head(5))
  cat("\n")

  return(list(
    global    = global,
    nodes     = node_df,
    community = comm
  ))
}


# =============================================================================
# SECTION 4 — NETWORK VISUALISATION
# =============================================================================

#' Visualise the co-occurrence network with multiple layout options.
#'
#' @param g           An igraph object.
#' @param node_df     Data frame from calculate_network_topology().
#' @param layout      Layout algorithm: "fr", "kk", "dh", "stress". Default = "fr".
#' @param colour_by   Node attribute for colouring: "module", "phylum", "degree",
#'                    "hub_score". Default = "module".
#' @param size_by     Node attribute for sizing: "degree", "hub_score",
#'                    "betweenness". Default = "degree".
#' @param show_labels Whether to show node labels. Default = TRUE.
#' @param label_hubs  Whether to label only hub taxa. Default = TRUE.
#' @param edge_alpha  Edge transparency. Default = 0.4.
#' @return A ggraph plot.

plot_network <- function(g,
                          node_df,
                          layout      = "fr",
                          colour_by   = "module",
                          size_by     = "degree",
                          show_labels = TRUE,
                          label_hubs  = TRUE,
                          edge_alpha  = 0.4) {

  cat("=== Network visualisation ===\n")

  # Attach node attributes to igraph
  for (col in colnames(node_df)) {
    if (col != "taxon" && col %in% names(vertex_attr(g)) == FALSE) {
      g <- set_vertex_attr(g, col, value = node_df[[col]][
        match(V(g)$name, node_df$taxon)
      ])
    }
  }

  # Convert to tidygraph
  tg <- as_tbl_graph(g)

  # Layout
  set.seed(42)
  lay <- switch(layout,
                fr     = create_layout(tg, layout = "fr"),
                kk     = create_layout(tg, layout = "kk"),
                dh     = create_layout(tg, layout = "dh"),
                stress = create_layout(tg, layout = "stress"),
                create_layout(tg, layout = "fr"))

  # Determine colour palette
  colour_vals <- node_df[[colour_by]][match(V(g)$name, node_df$taxon)]

  if (colour_by == "module") {
    n_mods  <- n_distinct(colour_vals, na.rm = TRUE)
    pal     <- if (n_mods <= 8) brewer.pal(max(3, n_mods), "Set2") else
               colorRampPalette(brewer.pal(8, "Set2"))(n_mods)
  } else if (colour_by == "phylum" && "Phylum" %in% colnames(node_df)) {
    phyla   <- unique(na.omit(node_df$Phylum))
    n_phyla <- length(phyla)
    pal     <- if (n_phyla <= 12) brewer.pal(max(3, n_phyla), "Paired") else
               colorRampPalette(brewer.pal(12, "Paired"))(n_phyla)
    names(pal) <- phyla
  }

  p <- ggraph(lay) +
    geom_edge_link(
      aes(edge_colour = color_e,
          edge_width  = weight),
      alpha = edge_alpha,
      show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    scale_edge_colour_identity()

  # Nodes
  if (colour_by %in% c("degree", "hub_score", "betweenness", "eigenvector")) {
    p <- p +
      geom_node_point(
        aes(size   = .data[[size_by]],
            colour = .data[[colour_by]]),
        alpha = 0.9
      ) +
      scale_colour_distiller(palette = "YlOrRd", direction = 1,
                              name = str_to_title(colour_by)) +
      scale_size(range = c(2, 10), name = str_to_title(size_by))
  } else {
    p <- p +
      geom_node_point(
        aes(size   = .data[[size_by]],
            fill   = as.factor(.data[[colour_by]])),
        shape = 21, colour = "white", stroke = 0.4, alpha = 0.9
      ) +
      scale_fill_brewer(palette = "Set2",
                        name   = str_to_title(colour_by)) +
      scale_size(range = c(2, 10), name = str_to_title(size_by))
  }

  # Labels
  if (show_labels) {
    label_data <- if (label_hubs) {
      lay[node_df$is_hub[match(lay$name, node_df$taxon)] %in% TRUE, ]
    } else lay

    if (nrow(label_data) > 0) {
      p <- p +
        geom_node_text(
          data      = label_data,
          aes(label = name),
          size      = 2.8,
          fontface  = "italic",
          repel     = TRUE,
          max.overlaps = 20,
          bg.color  = "white",
          bg.r      = 0.1
        )
    }
  }

  p <- p +
    labs(
      title    = paste0("Co-occurrence network - ", layout, " layout"),
      subtitle = paste0(vcount(g), " nodes | ", ecount(g), " edges | ",
                        "Coloured by ", colour_by),
      caption  = "Green edges = positive correlation | Red edges = negative correlation"
    ) +
    theme_network()

  cat("  Nodes:", vcount(g), "| Edges:", ecount(g), "\n\n")
  return(p)
}


#' Plot the network with nodes sized by abundance and coloured by phylum.
#' A cleaner alternative layout for publications.

plot_network_phylum <- function(g, node_df, ps = NULL) {

  # Add phylum info
  if (!is.null(ps) && !"Phylum" %in% colnames(node_df)) {
    tax_df <- tryCatch({
      as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE) %>%
        rownames_to_column("taxon") %>% select(taxon, Phylum)
    }, error = function(e) NULL)
    if (!is.null(tax_df)) node_df <- left_join(node_df, tax_df, by = "taxon")
  }

  node_df$Phylum[is.na(node_df$Phylum)] <- "Unknown"

  # Keep top 8 phyla, merge rest to "Other"
  top_phyla <- node_df %>%
    dplyr::count(Phylum, sort = TRUE) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::pull(Phylum)
  node_df$phylum_plot <- ifelse(node_df$Phylum %in% top_phyla,
                                 node_df$Phylum, "Other")

  for (col in c("degree", "hub_score", "phylum_plot")) {
    g <- set_vertex_attr(g, col, value = node_df[[col]][
      match(V(g)$name, node_df$taxon)
    ])
  }

  tg  <- as_tbl_graph(g)
  set.seed(42)
  lay <- create_layout(tg, layout = "fr")

  phyla_all <- unique(node_df$phylum_plot)
  n_phyla   <- length(phyla_all)
  pal       <- c(brewer.pal(min(8, n_phyla), "Set2"),
                 if (n_phyla > 8) rep("grey70", n_phyla - 8))
  names(pal) <- phyla_all

  p <- ggraph(lay) +
    geom_edge_link(
      aes(edge_colour = color_e, edge_width = weight),
      alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.2, 1.2), guide = "none") +
    scale_edge_colour_identity() +
    geom_node_point(
      aes(size = degree, fill = phylum_plot),
      shape = 21, colour = "white", stroke = 0.3, alpha = 0.9
    ) +
    scale_fill_manual(values = pal, name = "Phylum") +
    scale_size(range = c(2, 12), name = "Degree") +
    geom_node_text(
      data = lay[node_df$is_hub[match(lay$name, node_df$taxon)] %in% TRUE, ],
      aes(label = name), size = 2.8, fontface = "italic",
      repel = TRUE, max.overlaps = 20, bg.color = "white", bg.r = 0.1
    ) +
    labs(
      title    = "Co-occurrence network - phylum composition",
      subtitle = paste0(vcount(g), " nodes | ", ecount(g), " edges"),
      caption  = "Node size = degree | Hub taxa labelled"
    ) +
    theme_network()

  return(p)
}


# =============================================================================
# SECTION 5 — HUB TAXA ANALYSIS
# =============================================================================

#' Characterise hub taxa and visualise their connectivity patterns.
#'
#' @param g       An igraph object.
#' @param node_df Data frame from calculate_network_topology().
#' @param top_n   Number of hub taxa to show. Default = 20.
#' @return A list: hub table and plots.

analyse_hub_taxa <- function(g, node_df, top_n = 20) {

  cat("=== Hub taxa analysis ===\n")
  node_df <- node_df %>% mutate(across(where(is.list), ~as.vector(unlist(.))))

  hub_df <- node_df %>%
    arrange(desc(hub_score)) %>%
    slice_head(n = top_n) %>%
    mutate(taxon = fct_reorder(taxon, hub_score))

  cat("  Top", nrow(hub_df), "hub taxa:\n")
  print(hub_df %>% select(taxon, degree, betweenness, eigenvector, hub_score,
                            any_of("Phylum")) %>% head(10))
  cat("\n")

  hub_df$fill_var <- if ("Phylum" %in% colnames(hub_df)) hub_df$Phylum else "Hub"

  # --- Plot 1: Hub score ranking -------------------------------------------
  p_rank <- ggplot(hub_df,
                   aes(x = hub_score, y = taxon,
                       fill = if ("Phylum" %in% colnames(hub_df))
                                Phylum else "Hub")) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    scale_fill_brewer(palette = "Set2",
                      name = "Phylum") +
    labs(
      title    = paste0("Top ", top_n, " hub taxa"),
      subtitle = "Composite hub score (degree + betweenness + eigenvector centrality)",
      x        = "Hub score (z-scored composite)",
      y        = NULL
    ) +
    theme_microbiome() +
    theme(axis.text.y = element_text(face = "italic", size = 8))

  # --- Plot 2: Centrality scatter ------------------------------------------
  p_scatter <- ggplot(node_df,
                      aes(x = degree, y = betweenness,
                          colour = is_hub, size = eigenvector)) +
    geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(
      data     = node_df %>% filter(is_hub),
      aes(label = taxon),
      size     = 2.5, fontface = "italic",
      max.overlaps = 15, segment.size = 0.3
    ) +
    scale_colour_manual(
      values = c(`TRUE` = "#e74c3c", `FALSE` = "#3498db"),
      labels = c(`TRUE` = "Hub", `FALSE` = "Non-hub"),
      name   = "Status"
    ) +
    scale_size(range = c(1.5, 7), name = "Eigenvector") +
    labs(
      title    = "Node centrality - degree vs betweenness",
      subtitle = "Hub taxa highlighted | Size = eigenvector centrality",
      x        = "Degree (number of connections)",
      y        = "Betweenness centrality (normalised)"
    ) +
    theme_microbiome()

  # --- Plot 3: Degree distribution (power law check) -----------------------
  deg_dist <- data.frame(
    degree = as.numeric(igraph::degree(g))
  ) %>%
    dplyr::count(degree) %>%
    dplyr::mutate(freq = n / sum(n))

  p_degree <- ggplot(deg_dist, aes(x = degree, y = freq)) +
    geom_col(fill = "#3498db", alpha = 0.75, width = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "#e74c3c",
                linetype = "dashed", linewidth = 0.8) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title    = "Degree distribution",
      subtitle = "Red line = linear fit (log-log for power law check)",
      x        = "Degree (k)",
      y        = "Frequency P(k)"
    ) +
    theme_microbiome()

  combined <- (p_rank | p_scatter) / p_degree +
    plot_annotation(title = "Hub Taxa Analysis")

  return(list(
    hub_df   = hub_df,
    plot     = combined,
    p_rank   = p_rank,
    p_scatter = p_scatter,
    p_degree = p_degree
  ))
}


# =============================================================================
# SECTION 6 — COMMUNITY DETECTION AND MODULE ANALYSIS
# =============================================================================

#' Characterise community modules identified by Louvain algorithm.
#'
#' @param g         An igraph object.
#' @param node_df   Data frame from calculate_network_topology().
#' @param ps        A phyloseq object (for taxonomic composition of modules).
#' @param rank      Taxonomic rank for composition. Default = "Phylum".
#' @return A list: module summary, composition plot, and module network.

analyse_modules <- function(g, node_df, ps = NULL, rank = "Phylum") {

  cat("=== Module (community) analysis ===\n")
  node_df <- node_df %>% mutate(across(where(is.list), ~as.vector(unlist(.))))

  modules <- node_df %>%
    dplyr::count(module, name = "n_taxa") %>%
    dplyr::arrange(dplyr::desc(n_taxa))

  cat("  Modules detected:", nrow(modules), "\n")
  cat("  Module sizes:", paste(modules$n_taxa, collapse = ", "), "\n\n")

  # Per-module density and within-module hub
  mod_stats <- lapply(unique(node_df$module), function(m) {
    taxa_in_mod <- node_df$taxon[node_df$module == m]
    g_sub       <- induced_subgraph(g, vids = taxa_in_mod)
    top_hub     <- node_df %>%
      filter(module == m) %>%
      arrange(desc(hub_score)) %>%
      slice_head(n = 1) %>%
      pull(taxon)

    data.frame(
      module          = m,
      n_taxa          = vcount(g_sub),
      internal_edges  = ecount(g_sub),
      internal_density = round(edge_density(g_sub), 4),
      top_hub         = top_hub,
      stringsAsFactors = FALSE
    )
  })
  mod_stats_df <- bind_rows(mod_stats) %>% arrange(desc(n_taxa))
  print(mod_stats_df)
  cat("\n")

  # --- Phylum composition per module ----------------------------------------
  if (!is.null(ps) && rank %in% tryCatch(colnames(as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)), error = function(e) "")) {
    if (rank %in% colnames(node_df)) {
      node_comp <- node_df
    } else {
      tax_df <- tryCatch(as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column("taxon") %>%
        dplyr::select(taxon, dplyr::all_of(rank)), error = function(e) NULL)
      
      node_comp <- dplyr::left_join(node_df, tax_df, by = "taxon")
    }
    
    node_comp[[rank]][is.na(node_comp[[rank]])] <- "Unknown"

    top_phyla <- node_comp %>%
      dplyr::count(.data[[rank]], sort = TRUE) %>%
      dplyr::slice_head(n = 8) %>%
      dplyr::pull(.data[[rank]])
    node_comp[[paste0(rank, "_plot")]] <- ifelse(
      node_comp[[rank]] %in% top_phyla, node_comp[[rank]], "Other"
    )

    comp_df <- node_comp %>%
      dplyr::count(module, .data[[paste0(rank, "_plot")]]) %>%
      group_by(module) %>%
      mutate(pct = n / sum(n) * 100) %>%
      ungroup() %>%
      mutate(module = paste0("Module ", module))

    n_phyla <- n_distinct(comp_df[[paste0(rank, "_plot")]])
    pal     <- c(brewer.pal(min(8, n_phyla), "Paired"),
                 rep("grey70", max(0, n_phyla - 8)))

    p_comp <- ggplot(comp_df,
                     aes(x = module, y = pct,
                         fill = .data[[paste0(rank, "_plot")]])) +
      geom_col(position = "stack", width = 0.7, colour = "white",
               linewidth = 0.3) +
      scale_fill_manual(values = pal, name = rank) +
      scale_y_continuous(labels = label_number(suffix = "%"),
                         expand = c(0, 0)) +
      labs(
        title    = paste0(rank, " composition per network module"),
        subtitle = paste0(nrow(mod_stats_df), " modules detected"),
        x        = "Module",
        y        = "Percentage of taxa"
      ) +
      theme_microbiome() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    p_comp <- NULL
  }

  return(list(
    module_stats = mod_stats_df,
    plot_comp    = p_comp
  ))
}


# =============================================================================
# SECTION 7 — DIFFERENTIAL NETWORK ANALYSIS
# =============================================================================

#' Compare two networks (e.g. healthy vs disease) and identify
#' differential edges and connectivity changes.
#'
#' @param ps          A phyloseq object.
#' @param group_var   Metadata variable defining the two groups.
#' @param group1      Name of first group.
#' @param group2      Name of second group.
#' @param method      Correlation method. Default = "spearman".
#' @param cor_threshold Correlation threshold. Default = 0.6.
#' @param p_threshold  P-value threshold. Default = 0.05.
#' @return A list: both networks, differential analysis, and comparison plots.

compare_networks <- function(ps,
                              group_var     = "group",
                              group1        = NULL,
                              group2        = NULL,
                              method        = "spearman",
                              cor_threshold = 0.6,
                              p_threshold   = 0.05) {

  cat("=== Differential network analysis ===\n")
  cat("  Group 1:", group1, "| Group 2:", group2, "\n")

  groups <- unique(as.character(sample_data(ps)[[group_var]]))
  if (is.null(group1)) group1 <- groups[1]
  if (is.null(group2)) group2 <- groups[2]

  # Build separate networks
  build_group_network <- function(grp) {
    meta_df <- data.frame(phyloseq::sample_data(ps))
    
    keep_samples <- rownames(meta_df)[
      as.character(meta_df[[group_var]]) == grp
    ]
    
    ps_sub <- phyloseq::prune_samples(keep_samples, ps_sub <- ps)
    
    # Drop taxa that are zero in this group
    ps_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_sub) > 0, ps_sub)
    
    cat("    Samples:", phyloseq::nsamples(ps_sub),
        "| Taxa:", phyloseq::ntaxa(ps_sub), "\n")
    
    net <- build_network(
      ps_sub,
      method = method,
      cor_threshold = cor_threshold,
      p_threshold = p_threshold
    )
    
    topo <- calculate_network_topology(net$graph)
    
    list(network = net, topology = topo, group = grp)
  }

  cat("  Building", group1, "network...\n")
  net1 <- build_group_network(group1)
  cat("  Building", group2, "network...\n")
  net2 <- build_group_network(group2)

  # --- Compare global topology ----------------------------------------------
  comparison_df <- data.frame(
    Metric  = c("Nodes", "Edges", "Density", "Mean degree",
                "Clustering coefficient", "Avg path length",
                "Modularity", "N modules", "% positive edges"),
    Group1  = c(net1$topology$global$n_nodes,
                net1$topology$global$n_edges,
                net1$topology$global$density,
                net1$topology$global$mean_degree,
                net1$topology$global$clustering_coef,
                net1$topology$global$avg_path_length,
                net1$topology$global$modularity,
                net1$topology$global$n_modules,
                net1$topology$global$pct_positive),
    Group2  = c(net2$topology$global$n_nodes,
                net2$topology$global$n_edges,
                net2$topology$global$density,
                net2$topology$global$mean_degree,
                net2$topology$global$clustering_coef,
                net2$topology$global$avg_path_length,
                net2$topology$global$modularity,
                net2$topology$global$n_modules,
                net2$topology$global$pct_positive),
    stringsAsFactors = FALSE
  )
  colnames(comparison_df)[2:3] <- c(group1, group2)

  cat("\n  Global topology comparison:\n")
  print(comparison_df)
  cat("\n")

  # --- Differential edges ---------------------------------------------------
  # Edges present in one network but not the other
  edges1 <- as_edgelist(net1$network$graph) %>%
    as.data.frame() %>% setNames(c("from", "to")) %>%
    mutate(edge_id = paste(pmin(from, to), pmax(from, to), sep = "--"),
           present_in = group1)
  edges2 <- as_edgelist(net2$network$graph) %>%
    as.data.frame() %>% setNames(c("from", "to")) %>%
    mutate(edge_id = paste(pmin(from, to), pmax(from, to), sep = "--"),
           present_in = group2)

  unique_to_g1 <- edges1 %>% filter(!edge_id %in% edges2$edge_id)
  unique_to_g2 <- edges2 %>% filter(!edge_id %in% edges1$edge_id)
  shared_edges <- edges1 %>% filter(edge_id %in% edges2$edge_id)

  cat("  Edges unique to", group1, ":", nrow(unique_to_g1), "\n")
  cat("  Edges unique to", group2, ":", nrow(unique_to_g2), "\n")
  cat("  Shared edges:", nrow(shared_edges), "\n\n")

  # --- Comparison bar plot -------------------------------------------------
  comp_long <- comparison_df %>%
    pivot_longer(-Metric, names_to = "Group", values_to = "Value") %>%
    mutate(Value = as.numeric(Value))

  p_compare <- ggplot(comp_long,
                      aes(x = Group, y = Value, fill = Group)) +
    geom_col(width = 0.6, alpha = 0.85) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title    = "Network topology comparison",
      subtitle = paste0(group1, " vs ", group2),
      x        = NULL, y = "Value"
    ) +
    theme_microbiome() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- Edge summary Venn-like bar -------------------------------------------
  edge_summary <- data.frame(
    category = c(paste0("Unique to ", group1),
                 paste0("Unique to ", group2),
                 "Shared"),
    count    = c(nrow(unique_to_g1), nrow(unique_to_g2),
                 nrow(shared_edges))
  )

  p_edges <- ggplot(edge_summary, aes(x = category, y = count,
                                       fill = category)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = count), vjust = -0.4, fontface = "bold", size = 4) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Edge distribution between networks",
      x = NULL, y = "Number of edges"
    ) +
    theme_microbiome()

  return(list(
    net1          = net1,
    net2          = net2,
    comparison    = comparison_df,
    unique_g1     = unique_to_g1,
    unique_g2     = unique_to_g2,
    shared        = shared_edges,
    p_compare     = p_compare,
    p_edges       = p_edges
  ))
}


# =============================================================================
# SECTION 8 — RANDOM NETWORK COMPARISON
# =============================================================================

#' Compare observed network properties against random Erdos-Renyi networks.
#'
#' Tests whether the observed network has small-world properties
#' (high clustering, short path lengths) relative to random expectation.
#'
#' @param g           An igraph object.
#' @param n_random    Number of random networks to generate. Default = 100.
#' @param seed        Random seed. Default = 42.
#' @return A list: comparison table and plot.

compare_random_network <- function(g, n_random = 100, seed = 42) {

  cat("=== Random network comparison ===\n")
  cat("  Generating", n_random, "random networks...\n")

  n_nodes <- vcount(g)
  n_edges <- ecount(g)

  # Observed metrics
  obs <- list(
    clustering  = igraph::transitivity(g, type = "global"),
    path_length = igraph::mean_distance(g, directed = FALSE),
    degree_sd   = sd(igraph::degree(g))
  )

  # Random networks with same n_nodes and n_edges
  set.seed(seed)
  rand_metrics <- lapply(seq_len(n_random), function(i) {
    g_rand <- erdos.renyi.game(n_nodes, n_edges, type = "gnm")
    if (!is_connected(g_rand)) {
      # Keep only largest connected component
      comp   <- components(g_rand)
      g_rand <- induced_subgraph(g_rand,
                                  which(comp$membership == which.max(comp$csize)))
    }
    list(
      clustering  = transitivity(g_rand, type = "global"),
      path_length = mean_distance(g_rand, directed = FALSE)
    )
  })

  rand_df <- data.frame(
    clustering  = sapply(rand_metrics, `[[`, "clustering"),
    path_length = sapply(rand_metrics, `[[`, "path_length")
  )

  # Small-world index: C_obs/C_rand >> 1 and L_obs/L_rand ≈ 1
  sw_index <- (obs$clustering / mean(rand_df$clustering)) /
              (obs$path_length / mean(rand_df$path_length))

  cat("  Observed clustering:", round(obs$clustering, 4),
      "| Random mean:", round(mean(rand_df$clustering), 4), "\n")
  cat("  Observed path length:", round(obs$path_length, 4),
      "| Random mean:", round(mean(rand_df$path_length), 4), "\n")
  cat("  Small-world index (Sigma):", round(sw_index, 3))
  if (sw_index > 1) cat(" OK: Small-world network detected\n\n")
  else cat(" - Not small-world\n\n")

  # Plot distributions with observed value
  rand_long <- rand_df %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    mutate(metric = ifelse(metric == "clustering", "Clustering coefficient",
                           "Average path length"))

  obs_df <- data.frame(
    metric = c("Clustering coefficient", "Average path length"),
    value  = c(obs$clustering, obs$path_length)
  )

  p <- ggplot(rand_long, aes(x = value)) +
    geom_histogram(bins = 25, fill = "#3498db", colour = "white",
                   alpha = 0.75) +
    geom_vline(data = obs_df, aes(xintercept = value),
               colour = "#e74c3c", linewidth = 1.2, linetype = "dashed") +
    geom_text(data = obs_df,
              aes(x = value, y = Inf, label = paste0("Observed: ", round(value, 4))),
              vjust = 1.5, hjust = -0.1, colour = "#e74c3c", size = 3.5) +
    facet_wrap(~ metric, scales = "free", nrow = 1) +
    labs(
      title    = paste0("Observed vs random network (n = ", n_random, ")"),
      subtitle = paste0("Small-world index Sigma = ", round(sw_index, 3),
                        if (sw_index > 1) " OK: Small-world" else ""),
      x        = "Value",
      y        = "Count"
    ) +
    theme_microbiome()

  return(list(
    observed    = obs,
    random      = rand_df,
    sw_index    = sw_index,
    plot        = p
  ))
}


# =============================================================================
# SECTION 9 — COMPLETE NETWORK ANALYSIS WORKFLOW
# =============================================================================

#' Run the complete network analysis pipeline.
#'
#' @param ps          A filtered phyloseq object (raw counts).
#' @param group_var   Metadata variable for differential network analysis.
#' @param rank        Taxonomic rank. Default = "Genus".
#' @param method      Network method. Default = "spearman".
#' @param cor_threshold Correlation threshold. Default = 0.6.
#' @param output_dir  Directory to save all outputs.
#' @return A named list of all results.

run_network_analysis <- function(ps,
                                  group_var     = NULL,
                                  rank          = "Genus",
                                  method        = "spearman",
                                  cor_threshold = 0.6,
                                  output_dir    = "network_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  NETWORK ANALYSIS PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  results <- list()

  # --- Prepare data ---------------------------------------------------------
  ps_net <- prepare_network_data(ps, rank = rank,
                                  min_prevalence = 0.20, max_taxa = 100)

  # --- Build network --------------------------------------------------------
  cat("--- Building full network ---\n")
  net_res <- build_network(ps_net, method = method,
                            cor_threshold = cor_threshold)
  results$network <- net_res

  # --- Topology -------------------------------------------------------------
  cat("--- Network topology ---\n")
  topo <- calculate_network_topology(net_res$graph, ps = ps_net, rank = rank)
  results$topology <- topo
  write.csv(topo$nodes,
            file.path(output_dir, "network_node_metrics.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(topo$global),
            file.path(output_dir, "network_global_topology.csv"),
            row.names = FALSE)

  # --- Network plots --------------------------------------------------------
  cat("--- Plot 1: Network (module colouring) ---\n")
  results$p_network_module <- plot_network(
    net_res$graph, topo$nodes,
    layout = "fr", colour_by = "module",
    size_by = "degree", label_hubs = TRUE
  )
  ggsave(file.path(output_dir, "01_network_modules.pdf"),
         results$p_network_module, width = 12, height = 10)

  cat("--- Plot 2: Network (phylum colouring) ---\n")
  results$p_network_phylum <- plot_network_phylum(
    net_res$graph, topo$nodes, ps = ps_net
  )
  ggsave(file.path(output_dir, "02_network_phylum.pdf"),
         results$p_network_phylum, width = 12, height = 10)

  # --- Hub analysis ---------------------------------------------------------
  cat("--- Plot 3: Hub taxa ---\n")
  hub_res <- analyse_hub_taxa(net_res$graph, topo$nodes, top_n = 20)
  results$p_hubs <- hub_res$plot
  write.csv(hub_res$hub_df,
            file.path(output_dir, "hub_taxa.csv"),
            row.names = FALSE)
  ggsave(file.path(output_dir, "03_hub_taxa.pdf"),
         hub_res$plot, width = 14, height = 10)

  # --- Module analysis ------------------------------------------------------
  cat("--- Plot 4: Module composition ---\n")
  mod_res <- analyse_modules(net_res$graph, topo$nodes,
                              ps = ps_net, rank = "Phylum")
  if (!is.null(mod_res$plot_comp)) {
    results$p_modules <- mod_res$plot_comp
    ggsave(file.path(output_dir, "04_module_composition.pdf"),
           mod_res$plot_comp, width = 12, height = 7)
  }
  write.csv(mod_res$module_stats,
            file.path(output_dir, "module_statistics.csv"),
            row.names = FALSE)

  # --- Random network comparison --------------------------------------------
  cat("--- Plot 5: Random network comparison ---\n")
  rand_res <- compare_random_network(net_res$graph, n_random = 100)
  results$p_random  <- rand_res$plot
  results$sw_index  <- rand_res$sw_index
  ggsave(file.path(output_dir, "05_random_network_comparison.pdf"),
         rand_res$plot, width = 12, height = 6)

  # --- Differential network (optional) -------------------------------------
  if (!is.null(group_var)) {
    groups <- unique(as.character(sample_data(ps)[[group_var]]))
    if (length(groups) >= 2) {
      cat("--- Plot 6: Differential network ---\n")
      diff_net <- compare_networks(
        ps_net, group_var = group_var,
        group1 = groups[1], group2 = groups[2],
        method = method, cor_threshold = cor_threshold
      )
      results$differential <- diff_net
      write.csv(diff_net$comparison,
                file.path(output_dir, "network_topology_comparison.csv"),
                row.names = FALSE)
      ggsave(file.path(output_dir, "06_topology_comparison.pdf"),
             diff_net$p_compare, width = 14, height = 10)
      ggsave(file.path(output_dir, "07_edge_comparison.pdf"),
             diff_net$p_edges, width = 8, height = 6)
    }
  }

  cat("\n", strrep("=", 60), "\n")
  cat("  NETWORK ANALYSIS PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  up to 7 PDF files\n")
  cat("  Tables: node metrics, global topology, hubs, modules\n")
  cat("  Small-world index:", round(rand_res$sw_index, 3), "\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_network_analysis(
#   ps            = ps,
#   group_var     = "disease_status",
#   rank          = "Genus",
#   method        = "spearman",
#   cor_threshold = 0.6,
#   output_dir    = "results/network"
# )

# --- Option B: Step by step --------------------------------------------------
# ps_net    <- prepare_network_data(ps, rank = "Genus", min_prevalence = 0.20)
# net_res   <- build_network(ps_net, method = "spearman", cor_threshold = 0.6)
# topo      <- calculate_network_topology(net_res$graph, ps = ps_net)
# p_net     <- plot_network(net_res$graph, topo$nodes, colour_by = "module")
# hub_res   <- analyse_hub_taxa(net_res$graph, topo$nodes)
# hub_res$hub_df %>% select(taxon, degree, hub_score) %>% head(10)

# --- Option C: Differential network only ------------------------------------
# diff_res  <- compare_networks(ps_net, group_var = "disease_status",
#                               group1 = "Healthy", group2 = "IBD")
# diff_res$comparison
# diff_res$p_compare
