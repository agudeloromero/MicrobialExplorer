
simulate_microbiome_data <- function(
    output_dir = "data/simulated",
    dataset_type = c("network_enriched", "basic"),
    seed = 42,
    n_samples = 60,
    n_taxa = 200,
    n_subjects = 20,
    n_timepoints = 3,
    add_tree = TRUE
) {
  dataset_type <- match.arg(dataset_type)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  set.seed(seed)
  
  # OTU table
  otu_sim <- matrix(0, nrow = n_taxa, ncol = n_samples)
  
  for (i in seq_len(round(n_taxa * 0.8))) {
    prevalence <- runif(1, 0.6, 1.0)
    present <- sample(n_samples, round(n_samples * prevalence))
    otu_sim[i, present] <- rnbinom(length(present), mu = 80, size = 0.5)
  }
  
  for (i in (round(n_taxa * 0.8) + 1):n_taxa) {
    prevalence <- runif(1, 0.01, 0.08)
    present <- sample(n_samples, max(1, round(n_samples * prevalence)))
    otu_sim[i, present] <- rnbinom(length(present), mu = 5, size = 0.3)
  }
  
  if (dataset_type == "network_enriched") {
    base_pattern <- rnbinom(n_samples, mu = 100, size = 0.5)
    for (i in 1:10) {
      noise <- rnbinom(n_samples, mu = 20, size = 0.5)
      otu_sim[i, ] <- round(base_pattern * runif(1, 0.5, 1.5) + noise)
    }
    
    base_pattern2 <- rnbinom(n_samples, mu = 80, size = 0.4)
    for (i in 11:20) {
      noise <- rnbinom(n_samples, mu = 15, size = 0.4)
      otu_sim[i, ] <- round(base_pattern2 * runif(1, 0.5, 1.5) + noise)
    }
  }
  
  otu_sim[otu_sim < 0] <- 0
  
  rownames(otu_sim) <- paste0("Taxon_", stringr::str_pad(seq_len(n_taxa), 3, pad = "0"))
  colnames(otu_sim) <- paste0("Sample_", stringr::str_pad(seq_len(n_samples), 3, pad = "0"))
  
  # Metadata
  meta_sim <- data.frame(
    disease_status = rep(c("Healthy", "IBD"), each = n_samples / 2),
    age = round(rnorm(n_samples, 45, 12)),
    bmi = round(rnorm(n_samples, 25, 4), 1),
    crp = round(rexp(n_samples, 0.3), 2),
    calprotectin = round(rexp(n_samples, 0.01), 1),
    timepoint = rep(c(0, 4, 8), times = n_subjects),
    subject_id = rep(
      paste0("SUBJ_", stringr::str_pad(seq_len(n_subjects), 2, pad = "0")),
      each = n_timepoints
    ),
    sex = sample(c("M", "F"), n_samples, replace = TRUE),
    row.names = colnames(otu_sim),
    stringsAsFactors = FALSE
  )
  
  # Taxonomy
  phyla <- c(
    "Firmicutes", "Bacteroidota", "Proteobacteria",
    "Actinobacteria", "Verrucomicrobiota"
  )
  
  tax_sim <- data.frame(
    Kingdom = "Bacteria",
    Phylum = sample(
      phyla, n_taxa, replace = TRUE,
      prob = c(0.40, 0.30, 0.15, 0.10, 0.05)
    ),
    Class = paste0("Class_", seq_len(n_taxa)),
    Order = paste0("Order_", seq_len(n_taxa)),
    Family = paste0("Family_", seq_len(n_taxa)),
    Genus = paste0("Genus_", seq_len(n_taxa)),
    Species = paste0("Species_", seq_len(n_taxa)),
    row.names = rownames(otu_sim),
    stringsAsFactors = FALSE
  )
  
  ps_sim <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_sim, taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(tax_sim)),
    phyloseq::sample_data(meta_sim)
  )
  
  if (add_tree) {
    fake_tree <- ape::rtree(phyloseq::ntaxa(ps_sim))
    fake_tree$tip.label <- phyloseq::taxa_names(ps_sim)
    ps_sim <- phyloseq::merge_phyloseq(
      ps_sim,
      phyloseq::phy_tree(fake_tree)
    )
  }
  
  # Save files
  otu_file <- file.path(output_dir, paste0(dataset_type, "_otu.csv"))
  tax_file <- file.path(output_dir, paste0(dataset_type, "_taxonomy.csv"))
  meta_file <- file.path(output_dir, paste0(dataset_type, "_metadata.csv"))
  rds_file <- file.path(output_dir, paste0(dataset_type, "_phyloseq.rds"))
  tree_file <- file.path(output_dir, paste0(dataset_type, "_tree.nwk"))
  
  write.csv(otu_sim, otu_file)
  write.csv(tax_sim, tax_file)
  write.csv(meta_sim, meta_file)
  saveRDS(ps_sim, rds_file)
  
  if (add_tree) {
    ape::write.tree(phyloseq::phy_tree(ps_sim), tree_file)
  } else {
    tree_file <- NULL
  }
  
  cat("Simulated dataset created:\n")
  cat("  Type:", dataset_type, "\n")
  cat("  Samples:", phyloseq::nsamples(ps_sim), "\n")
  cat("  Taxa:", phyloseq::ntaxa(ps_sim), "\n")
  cat("  OTU:", otu_file, "\n")
  cat("  Taxonomy:", tax_file, "\n")
  cat("  Metadata:", meta_file, "\n")
  cat("  Tree:", tree_file, "\n")
  cat("  Phyloseq RDS:", rds_file, "\n")
  
  return(list(
    ps = ps_sim,
    otu = otu_sim,
    taxonomy = tax_sim,
    metadata = meta_sim,
    files = list(
      otu = otu_file,
      taxonomy = tax_file,
      metadata = meta_file,
      tree = tree_file,
      phyloseq = rds_file
    )
  ))
}
