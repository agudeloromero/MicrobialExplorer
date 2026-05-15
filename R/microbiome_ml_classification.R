# =============================================================================
# MICROBIOME MACHINE LEARNING AND CLASSIFICATION
# =============================================================================
# Description : Comprehensive ML pipeline for microbiome classification
#               covering random forest, cross-validation, feature importance,
#               ROC analysis, biomarker discovery, and model evaluation
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Trained models, ROC curves, feature importance plots,
#               confusion matrices, biomarker panels, SHAP-like explanations
# Author      : Patricia
# Dependencies: phyloseq, randomForest, caret, pROC, ggplot2, dplyr, tidyr,
#               patchwork, scales, tibble, RColorBrewer, vegan, e1071,
#               glmnet, xgboost (optional), DALEX (optional)
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(tibble)
  library(RColorBrewer)
  library(vegan)
  library(e1071)
  library(glmnet)
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


# =============================================================================
# SECTION 1 — DATA PREPARATION FOR ML
# =============================================================================

#' Prepare a phyloseq object for machine learning classification.
#'
#' Applies CLR transformation, handles class imbalance, performs
#' train/test splitting, and optionally applies feature pre-selection
#' to reduce dimensionality before modelling.
#'
#' @param ps              A phyloseq object (raw counts).
#' @param group_var       Target variable (class label) for classification.
#' @param rank            Taxonomic rank to agglomerate to. Default = "Genus".
#' @param transform       Feature transformation: "clr", "relative", or "log10".
#' @param test_fraction   Fraction of samples held out for testing. Default = 0.2.
#' @param min_prevalence  Minimum prevalence for feature inclusion. Default = 0.10.
#' @param max_features    Maximum features to keep. Default = 200.
#' @param balance_classes Whether to apply SMOTE-like upsampling. Default = FALSE.
#' @param seed            Random seed. Default = 42.
#' @return A list: training set, test set, feature names, class labels.

prepare_ml_data <- function(ps,
                             group_var       = "group",
                             rank            = "Genus",
                             transform       = "clr",
                             test_fraction   = 0.2,
                             min_prevalence  = 0.10,
                             max_features    = 200,
                             balance_classes = FALSE,
                             seed            = 42) {

  cat("=== Preparing data for machine learning ===\n")

  # --- Agglomerate and filter -----------------------------------------------
  ps_agg   <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  dup_idx  <- duplicated(new_names)
  new_names <- make.unique(new_names, sep = "_dup")
  taxa_names(ps_agg) <- new_names

  otu_mat <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Prevalence filter
  prev    <- rowSums(otu_mat > 0) / ncol(otu_mat)
  keep    <- prev >= min_prevalence
  otu_mat <- otu_mat[keep, ]

  # Trim to max_features by variance
  if (nrow(otu_mat) > max_features) {
    feat_var <- apply(otu_mat, 1, var)
    top_idx  <- order(feat_var, decreasing = TRUE)[seq_len(max_features)]
    otu_mat  <- otu_mat[top_idx, ]
  }

  cat("  Features after filtering:", nrow(otu_mat), "\n")

  # --- Transformation -------------------------------------------------------
  if (transform == "clr") {
    feat_mat <- t(apply(otu_mat + 0.5, 2, function(x) log(x) - mean(log(x))))
  } else if (transform == "relative") {
    feat_mat <- t(apply(otu_mat, 2, function(x) x / sum(x)))
  } else if (transform == "log10") {
    feat_mat <- t(log10(otu_mat + 1))
  } else {
    feat_mat <- t(otu_mat)
  }

  # Samples as rows, features as columns
  feat_df  <- as.data.frame(feat_mat)
  labels   <- as.factor(sample_data(ps_agg)[[group_var]])
  feat_df$label <- labels

  cat("  Samples:", nrow(feat_df), "\n")
  cat("  Class distribution:\n")
  print(table(feat_df$label))

  # --- Train / test split ---------------------------------------------------
  set.seed(seed)
  train_idx <- createDataPartition(feat_df$label,
                                    p     = 1 - test_fraction,
                                    list  = FALSE)[, 1]
  train_df  <- feat_df[train_idx, ]
  test_df   <- feat_df[-train_idx, ]

  cat("\n  Training samples:", nrow(train_df),
      "| Test samples:", nrow(test_df), "\n\n")

  # --- Class balancing (upsampling minority class) --------------------------
  if (balance_classes) {
    set.seed(seed)
    train_df <- upSample(
      x     = train_df %>% select(-label),
      y     = train_df$label,
      yname = "label"
    )
    cat("  After upsampling:\n")
    print(table(train_df$label))
    cat("\n")
  }

  return(list(
    train      = train_df,
    test       = test_df,
    features   = setdiff(colnames(feat_df), "label"),
    labels     = levels(labels),
    n_classes  = nlevels(labels),
    group_var  = group_var
  ))
}


# =============================================================================
# SECTION 2 — RANDOM FOREST
# =============================================================================

#' Train a Random Forest classifier with cross-validation and tuning.
#'
#' @param ml_data       Output from prepare_ml_data().
#' @param n_trees       Number of trees. Default = 500.
#' @param cv_folds      Number of cross-validation folds. Default = 5.
#' @param cv_repeats    Number of CV repeats. Default = 3.
#' @param tune_mtry     Whether to tune mtry parameter. Default = TRUE.
#' @param seed          Random seed. Default = 42.
#' @return A list: trained model, CV results, and performance metrics.

train_random_forest <- function(ml_data,
                                 n_trees    = 500,
                                 cv_folds   = 5,
                                 cv_repeats = 3,
                                 tune_mtry  = TRUE,
                                 seed       = 42) {

  cat("=== Training Random Forest ===\n")

  train_df  <- ml_data$train
  features  <- ml_data$features
  X_train   <- train_df[, features]
  y_train   <- train_df$label

  # --- Cross-validation setup -----------------------------------------------
  ctrl <- trainControl(
    method          = "repeatedcv",
    number          = cv_folds,
    repeats         = cv_repeats,
    classProbs      = TRUE,
    summaryFunction = if (ml_data$n_classes == 2) twoClassSummary
                      else multiClassSummary,
    savePredictions = "final",
    verboseIter     = FALSE
  )

  # Tune grid
  mtry_vals  <- if (tune_mtry) {
    floor(c(sqrt(length(features)),
            sqrt(length(features)) * 0.5,
            sqrt(length(features)) * 1.5,
            sqrt(length(features)) * 2))
  } else {
    floor(sqrt(length(features)))
  }
  tune_grid  <- expand.grid(mtry = unique(pmax(1, mtry_vals)))

  cat("  CV folds:", cv_folds, "| Repeats:", cv_repeats,
      "| mtry values:", paste(tune_grid$mtry, collapse = ", "), "\n")

  set.seed(seed)
  rf_model <- caret::train(
    x          = X_train,
    y          = y_train,
    method     = "rf",
    ntree      = n_trees,
    tuneGrid   = tune_grid,
    trControl  = ctrl,
    metric     = if (ml_data$n_classes == 2) "ROC" else "Accuracy",
    importance = TRUE
  )

  cat("  Best mtry:", rf_model$bestTune$mtry, "\n")
  cat("  CV", if (ml_data$n_classes == 2) "AUC" else "Accuracy", ":",
      round(max(rf_model$results[[if (ml_data$n_classes == 2) "ROC" else "Accuracy"]]),
            4), "\n\n")

  # --- Test set evaluation --------------------------------------------------
  test_df   <- ml_data$test
  X_test    <- test_df[, features]
  y_test    <- test_df$label

  pred_class <- predict(rf_model, X_test)
  pred_prob  <- predict(rf_model, X_test, type = "prob")

  cm <- confusionMatrix(pred_class, y_test)
  cat("  Test set performance:\n")
  cat("    Accuracy:", round(cm$overall["Accuracy"], 4), "\n")
  if (ml_data$n_classes == 2) {
    cat("    Sensitivity:", round(cm$byClass["Sensitivity"], 4), "\n")
    cat("    Specificity:", round(cm$byClass["Specificity"], 4), "\n")
  }
  cat("\n")

  return(list(
    model        = rf_model,
    predictions  = pred_class,
    probabilities = pred_prob,
    confusion    = cm,
    cv_results   = rf_model$results,
    best_mtry    = rf_model$bestTune$mtry
  ))
}


# =============================================================================
# SECTION 3 — LASSO LOGISTIC REGRESSION
# =============================================================================

#' Train a LASSO logistic regression for sparse biomarker selection.
#'
#' LASSO performs automatic feature selection by shrinking coefficients
#' of uninformative features to zero. Ideal for identifying minimal
#' biomarker panels.
#'
#' @param ml_data     Output from prepare_ml_data().
#' @param cv_folds    Number of cross-validation folds. Default = 10.
#' @param alpha       Elastic net mixing parameter (1 = LASSO, 0 = Ridge).
#' @param seed        Random seed. Default = 42.
#' @return A list: trained model, selected features, and performance.

train_lasso <- function(ml_data,
                         cv_folds = 10,
                         alpha    = 1,
                         seed     = 42) {

  cat("=== Training LASSO regression ===\n")

  train_df <- ml_data$train
  features <- ml_data$features
  X_train  <- as.matrix(train_df[, features])
  y_train  <- train_df$label

  # Binary: encode as 0/1
  if (ml_data$n_classes == 2) {
    y_bin    <- as.integer(y_train) - 1
    family   <- "binomial"
  } else {
    y_bin    <- y_train
    family   <- "multinomial"
  }

  set.seed(seed)
  cv_lasso <- cv.glmnet(
    x         = X_train,
    y         = y_bin,
    family    = family,
    alpha     = alpha,
    nfolds    = cv_folds,
    type.measure = "auc"
  )

  best_lambda <- cv_lasso$lambda.1se
  cat("  Best lambda (1se):", round(best_lambda, 6), "\n")
  cat("  Lambda min AUC:",
      round(max(cv_lasso$cvm, na.rm = TRUE), 4), "\n")

  # Selected features at best lambda
  coef_mat  <- coef(cv_lasso, s = "lambda.1se")

  if (ml_data$n_classes == 2) {
    coef_df <- as.data.frame(as.matrix(coef_mat)) %>%
      rownames_to_column("feature") %>%
      setNames(c("feature", "coefficient")) %>%
      filter(feature != "(Intercept)", coefficient != 0) %>%
      arrange(desc(abs(coefficient)))
  } else {
    # Multinomial: combine coefficients across classes
    coef_df <- lapply(names(coef_mat), function(cls) {
      df <- as.data.frame(as.matrix(coef_mat[[cls]])) %>%
        rownames_to_column("feature") %>%
        setNames(c("feature", "coefficient")) %>%
        filter(feature != "(Intercept)", coefficient != 0) %>%
        mutate(class = cls)
      df
    }) %>% bind_rows() %>% arrange(desc(abs(coefficient)))
  }

  n_selected <- n_distinct(coef_df$feature)
  cat("  Features selected:", n_selected, "\n\n")

  # Test set evaluation
  test_df   <- ml_data$test
  X_test    <- as.matrix(test_df[, features])
  y_test    <- test_df$label

  pred_prob <- predict(cv_lasso, X_test, s = "lambda.1se", type = "response")

  if (ml_data$n_classes == 2) {
    pred_class <- factor(ifelse(pred_prob > 0.5,
                                 ml_data$labels[2],
                                 ml_data$labels[1]),
                          levels = ml_data$labels)
    cm <- confusionMatrix(pred_class, y_test)
    cat("  LASSO test accuracy:", round(cm$overall["Accuracy"], 4), "\n\n")
  } else {
    cm <- NULL
  }

  return(list(
    model           = cv_lasso,
    selected_features = coef_df,
    n_selected      = n_selected,
    best_lambda     = best_lambda,
    confusion       = cm,
    cv_object       = cv_lasso
  ))
}


# =============================================================================
# SECTION 4 — FEATURE IMPORTANCE
# =============================================================================

#' Extract and visualise feature importance from Random Forest and LASSO.
#'
#' @param rf_result   Output from train_random_forest().
#' @param lasso_result Output from train_lasso().
#' @param top_n       Number of top features to display. Default = 30.
#' @param ps          A phyloseq object (for adding taxonomy labels).
#' @param rank        Taxonomic rank label. Default = "Genus".
#' @return A list: importance data frames and combined plot.

plot_feature_importance <- function(rf_result    = NULL,
                                     lasso_result = NULL,
                                     top_n        = 30,
                                     ps           = NULL,
                                     rank         = "Genus") {

  cat("=== Feature importance ===\n")

  plots <- list()
  imp_tables <- list()

  # --- Random Forest importance ---------------------------------------------
  if (!is.null(rf_result)) {
    rf_imp <- varImp(rf_result$model, scale = TRUE)$importance

    if (ncol(rf_imp) > 1) {
      # Multi-class: mean importance across classes
      rf_imp$Overall <- rowMeans(rf_imp)
    }

    rf_imp_df <- rf_imp %>%
      rownames_to_column("feature") %>%
      arrange(desc(Overall)) %>%
      slice_head(n = top_n) %>%
      mutate(feature = fct_reorder(feature, Overall))

    # Add phylum annotation if available
    if (!is.null(ps)) {
      tax_df <- as.data.frame(tax_table(ps)) %>%
        rownames_to_column("feature") %>%
        select(feature, any_of("Phylum"))
      rf_imp_df <- left_join(rf_imp_df, tax_df, by = "feature")
    }

    imp_tables$rf <- rf_imp_df

    p_rf <- ggplot(rf_imp_df,
                   aes(x = Overall, y = feature,
                       fill = if ("Phylum" %in% colnames(rf_imp_df) &&
                                  !all(is.na(rf_imp_df$Phylum)))
                                Phylum else "RF")) +
      geom_col(width = 0.7, alpha = 0.85) +
      scale_fill_brewer(palette = "Set2",
                        name = if ("Phylum" %in% colnames(rf_imp_df)) "Phylum"
                               else NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(
        title    = paste0("Random Forest — top ", top_n, " features"),
        subtitle = "Mean decrease in accuracy (scaled 0–100)",
        x        = "Importance",
        y        = NULL
      ) +
      theme_microbiome() +
      theme(axis.text.y = element_text(face = "italic", size = 8))

    plots$rf <- p_rf
    cat("  RF importance: top feature =",
        as.character(rf_imp_df$feature[1]), "\n")
  }

  # --- LASSO coefficients ---------------------------------------------------
  if (!is.null(lasso_result) && nrow(lasso_result$selected_features) > 0) {
    lasso_df <- lasso_result$selected_features %>%
      filter("feature" %in% colnames(.) | TRUE) %>%
      arrange(desc(abs(coefficient))) %>%
      slice_head(n = top_n) %>%
      mutate(
        feature   = fct_reorder(feature, abs(coefficient)),
        direction = ifelse(coefficient > 0, "Positive", "Negative")
      )

    imp_tables$lasso <- lasso_df

    colour_map <- c("Positive" = "#e74c3c", "Negative" = "#3498db")

    p_lasso <- ggplot(lasso_df,
                      aes(x = coefficient, y = feature, fill = direction)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_vline(xintercept = 0, linewidth = 0.7, colour = "grey30") +
      scale_fill_manual(values = colour_map, name = "Direction") +
      labs(
        title    = paste0("LASSO — selected features (n=",
                          lasso_result$n_selected, ")"),
        subtitle = "Non-zero coefficients at lambda.1se",
        x        = "LASSO coefficient",
        y        = NULL
      ) +
      theme_microbiome() +
      theme(axis.text.y = element_text(face = "italic", size = 8))

    plots$lasso <- p_lasso
    cat("  LASSO selected features:", lasso_result$n_selected, "\n")
  }

  # --- Consensus: features selected by both methods -------------------------
  if (!is.null(rf_result) && !is.null(lasso_result)) {
    rf_top     <- as.character(imp_tables$rf$feature[seq_len(min(50, nrow(imp_tables$rf)))])
    lasso_feats <- unique(as.character(lasso_result$selected_features$feature))
    consensus   <- intersect(rf_top, lasso_feats)

    cat("  Consensus features (RF top 50 ∩ LASSO):", length(consensus), "\n")

    if (length(consensus) > 0) {
      cons_df <- imp_tables$rf %>%
        filter(as.character(feature) %in% consensus) %>%
        mutate(feature = fct_reorder(as.character(feature), Overall))

      p_cons <- ggplot(cons_df,
                       aes(x = Overall, y = feature,
                           fill = if ("Phylum" %in% colnames(cons_df) &&
                                      !all(is.na(cons_df$Phylum)))
                                    Phylum else "Consensus")) +
        geom_col(width = 0.65, alpha = 0.85) +
        geom_text(aes(x = 0.5, label = feature),
                  hjust = 0, size = 2.8, fontface = "italic",
                  colour = "white") +
        scale_fill_brewer(palette = "Set1",
                          name = if ("Phylum" %in% colnames(cons_df)) "Phylum"
                                 else NULL) +
        labs(
          title    = paste0("Consensus biomarkers (n = ", length(consensus), ")"),
          subtitle = "Selected by both Random Forest and LASSO",
          x        = "RF importance",
          y        = NULL
        ) +
        theme_microbiome() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

      plots$consensus <- p_cons
    }
  }

  cat("\n")

  # Arrange plots
  if (length(plots) == 3) {
    combined <- (plots$rf | plots$lasso) / plots$consensus +
      plot_annotation(title = "Feature Importance Analysis",
                      theme = theme(plot.title = element_text(size = 14, face = "bold")))
  } else if (length(plots) == 2) {
    combined <- plots[[1]] | plots[[2]]
  } else {
    combined <- plots[[1]]
  }

  return(list(
    plot       = combined,
    plots      = plots,
    importance = imp_tables
  ))
}


# =============================================================================
# SECTION 5 — ROC CURVE ANALYSIS
# =============================================================================

#' Generate ROC curves for one or multiple classifiers.
#'
#' @param model_results Named list of model result objects (each with $probabilities).
#' @param ml_data       Output from prepare_ml_data().
#' @param positive_class The positive class label. Default = second class level.
#' @return A list: AUC values, ROC objects, and plot.

plot_roc_curves <- function(model_results,
                              ml_data,
                              positive_class = NULL) {

  cat("=== ROC curve analysis ===\n")

  y_test  <- ml_data$test$label

  if (is.null(positive_class)) positive_class <- ml_data$labels[2]

  roc_list <- list()
  auc_df   <- data.frame()

  colours  <- brewer.pal(max(3, length(model_results)), "Set1")[
    seq_len(length(model_results))
  ]
  names(colours) <- names(model_results)

  for (mod_name in names(model_results)) {
    res      <- model_results[[mod_name]]

    if (is.null(res$probabilities)) next

    # Handle matrix or vector probabilities
    if (is.matrix(res$probabilities) || is.data.frame(res$probabilities)) {
      if (positive_class %in% colnames(res$probabilities)) {
        prob_pos <- res$probabilities[, positive_class]
      } else {
        prob_pos <- res$probabilities[, 2]
      }
    } else {
      prob_pos <- as.numeric(res$probabilities)
    }

    roc_obj       <- roc(y_test, prob_pos,
                          levels    = c(setdiff(ml_data$labels, positive_class),
                                        positive_class),
                          direction = "<",
                          quiet     = TRUE)
    roc_list[[mod_name]] <- roc_obj

    auc_val <- round(auc(roc_obj), 4)
    ci_val  <- ci.auc(roc_obj, conf.level = 0.95)

    auc_df <- rbind(auc_df, data.frame(
      model   = mod_name,
      auc     = auc_val,
      ci_lo   = round(ci_val[1], 4),
      ci_hi   = round(ci_val[3], 4),
      stringsAsFactors = FALSE
    ))

    cat(" ", mod_name, "AUC:", auc_val,
        "(95% CI:", round(ci_val[1], 3), "–", round(ci_val[3], 3), ")\n")
  }

  # --- Build ROC plot -------------------------------------------------------
  roc_df <- lapply(names(roc_list), function(mod_name) {
    roc_obj <- roc_list[[mod_name]]
    auc_val <- auc_df$auc[auc_df$model == mod_name]
    data.frame(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      model       = paste0(mod_name, " (AUC=", auc_val, ")"),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  roc_colours <- setNames(colours, unique(roc_df$model))

  p_roc <- ggplot(roc_df,
                  aes(x = 1 - specificity, y = sensitivity,
                      colour = model)) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "grey70", linewidth = 0.8) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    scale_colour_brewer(palette = "Set1", name = "Model (AUC)") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0.01, 0)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0.01, 0)) +
    coord_equal() +
    labs(
      title    = "ROC curves",
      subtitle = paste0("Positive class: ", positive_class,
                        " | n test = ", nrow(ml_data$test)),
      x        = "1 − Specificity (False Positive Rate)",
      y        = "Sensitivity (True Positive Rate)"
    ) +
    theme_microbiome()

  # --- AUC comparison bar ---------------------------------------------------
  p_auc <- ggplot(auc_df,
                  aes(x = fct_reorder(model, auc), y = auc,
                      fill = model)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = ci_lo, ymax = ci_hi),
      width = 0.2, linewidth = 0.7
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey50", linewidth = 0.7) +
    geom_text(aes(label = auc), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0),
                       labels = label_number(accuracy = 0.01)) +
    labs(
      title    = "AUC comparison",
      subtitle = "Error bars = 95% CI",
      x        = NULL,
      y        = "AUC"
    ) +
    theme_microbiome() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  combined <- p_roc | p_auc
  cat("\n")

  return(list(
    roc_objects = roc_list,
    auc_table   = auc_df,
    p_roc       = p_roc,
    p_auc       = p_auc,
    plot        = combined
  ))
}


# =============================================================================
# SECTION 6 — CONFUSION MATRIX VISUALISATION
# =============================================================================

#' Visualise confusion matrix as a heatmap with performance metrics.
#'
#' @param confusion_matrix  A confusionMatrix object from caret.
#' @param model_name        Label for the plot title.
#' @return A ggplot heatmap.

plot_confusion_matrix <- function(confusion_matrix,
                                   model_name = "Model") {

  cat("=== Confusion matrix visualisation ===\n")

  cm_table <- as.data.frame(confusion_matrix$table)

  # Normalise by true class (row percentages)
  cm_norm <- cm_table %>%
    group_by(Reference) %>%
    mutate(pct = round(Freq / sum(Freq) * 100, 1)) %>%
    ungroup()

  n_classes <- n_distinct(cm_norm$Prediction)
  max_freq  <- max(cm_norm$Freq)

  p <- ggplot(cm_norm,
              aes(x = Reference, y = Prediction, fill = pct)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = paste0(Freq, "\n(", pct, "%)")),
              size = 4, fontface = "bold",
              colour = ifelse(cm_norm$pct > 60, "white", "black")) +
    scale_fill_gradient2(
      low      = "white",
      mid      = "#3498db",
      high     = "#1a5276",
      midpoint = 50,
      limits   = c(0, 100),
      name     = "Row %"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = paste0("Confusion matrix — ", model_name),
      subtitle = paste0(
        "Accuracy: ",
        round(confusion_matrix$overall["Accuracy"] * 100, 1), "% | ",
        "Kappa: ",
        round(confusion_matrix$overall["Kappa"], 3)
      ),
      x        = "True label",
      y        = "Predicted label"
    ) +
    theme_microbiome() +
    theme(panel.border = element_rect(colour = "grey80", fill = NA))

  # Performance metrics panel
  if (!is.null(confusion_matrix$byClass)) {
    if (is.matrix(confusion_matrix$byClass)) {
      # Multi-class
      metrics_df <- as.data.frame(confusion_matrix$byClass) %>%
        rownames_to_column("class") %>%
        select(class, Sensitivity, Specificity,
               `Pos Pred Value`, `Balanced Accuracy`) %>%
        rename(Precision = `Pos Pred Value`,
               `Bal. Acc.` = `Balanced Accuracy`) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3)))
    } else {
      # Binary
      metrics_df <- data.frame(
        Metric = c("Sensitivity", "Specificity", "Precision",
                   "F1", "Balanced Accuracy"),
        Value  = round(c(
          confusion_matrix$byClass["Sensitivity"],
          confusion_matrix$byClass["Specificity"],
          confusion_matrix$byClass["Pos Pred Value"],
          confusion_matrix$byClass["F1"],
          confusion_matrix$byClass["Balanced Accuracy"]
        ), 4)
      )

      p_metrics <- ggplot(metrics_df, aes(x = Metric, y = Value,
                                           fill = Metric)) +
        geom_col(width = 0.65, alpha = 0.85) +
        geom_text(aes(label = Value), vjust = -0.4, size = 3.5,
                  fontface = "bold") +
        geom_hline(yintercept = 0.5, linetype = "dashed",
                   colour = "grey50") +
        scale_fill_brewer(palette = "Set2", guide = "none") +
        scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
        labs(
          title = "Performance metrics",
          x     = NULL, y = "Value"
        ) +
        theme_microbiome() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))

      return(list(
        plot         = p | p_metrics,
        p_cm         = p,
        p_metrics    = p_metrics
      ))
    }
  }

  return(list(plot = p, p_cm = p))
}


# =============================================================================
# SECTION 7 — CROSS-VALIDATION LEARNING CURVES
# =============================================================================

#' Plot learning curves to assess overfitting and data sufficiency.
#'
#' @param ml_data     Output from prepare_ml_data().
#' @param fractions   Training size fractions. Default = seq(0.2, 1.0, by=0.1).
#' @param cv_folds    CV folds for each size. Default = 5.
#' @param n_trees     RF trees. Default = 200.
#' @param seed        Random seed. Default = 42.
#' @return A ggplot learning curve.

plot_learning_curve <- function(ml_data,
                                 fractions = seq(0.2, 1.0, by = 0.1),
                                 cv_folds  = 5,
                                 n_trees   = 200,
                                 seed      = 42) {

  cat("=== Learning curve analysis ===\n")

  train_df <- ml_data$train
  features <- ml_data$features
  X_all    <- train_df[, features]
  y_all    <- train_df$label

  curve_results <- lapply(fractions, function(frac) {
    set.seed(seed)
    subset_idx <- createDataPartition(y_all, p = frac, list = FALSE)[, 1]
    X_sub      <- X_all[subset_idx, ]
    y_sub      <- y_all[subset_idx]

    ctrl <- trainControl(
      method      = "cv",
      number      = min(cv_folds, length(subset_idx) - 1),
      classProbs  = TRUE,
      verboseIter = FALSE,
      summaryFunction = if (ml_data$n_classes == 2) twoClassSummary
                        else multiClassSummary
    )

    set.seed(seed)
    mod <- tryCatch({
      caret::train(x = X_sub, y = y_sub, method = "rf",
            ntree = n_trees,
            tuneGrid = data.frame(mtry = floor(sqrt(ncol(X_sub)))),
            trControl = ctrl,
            metric = if (ml_data$n_classes == 2) "ROC" else "Accuracy")
    }, error = function(e) NULL)

    if (is.null(mod)) return(NULL)

    metric_col <- if (ml_data$n_classes == 2) "ROC" else "Accuracy"
    cv_score   <- max(mod$results[[metric_col]])

    # Training score
    pred_train <- predict(mod, X_sub)
    train_acc  <- sum(pred_train == y_sub) / length(y_sub)

    data.frame(
      n_samples  = length(subset_idx),
      fraction   = frac,
      cv_score   = round(cv_score, 4),
      train_acc  = round(train_acc, 4),
      stringsAsFactors = FALSE
    )
  })

  curve_df <- bind_rows(curve_results)
  cat("  Learning curve computed for", nrow(curve_df), "training sizes\n\n")

  metric_label <- if (ml_data$n_classes == 2) "AUC (ROC)" else "Accuracy"

  p <- ggplot(curve_df, aes(x = n_samples)) +
    geom_line(aes(y = cv_score, colour = "CV score"),
              linewidth = 1.2) +
    geom_point(aes(y = cv_score, colour = "CV score"), size = 3) +
    geom_line(aes(y = train_acc, colour = "Training score"),
              linewidth = 1.2, linetype = "dashed") +
    geom_point(aes(y = train_acc, colour = "Training score"), size = 3) +
    geom_ribbon(
      aes(ymin = cv_score, ymax = train_acc),
      fill = "#f39c12", alpha = 0.1
    ) +
    scale_colour_manual(
      values = c("CV score" = "#3498db", "Training score" = "#e74c3c"),
      name   = "Score type"
    ) +
    scale_y_continuous(limits = c(
      min(0.4, min(curve_df$cv_score) - 0.05),
      1.02
    )) +
    scale_x_continuous(breaks = pretty_breaks(n = 6)) +
    labs(
      title    = "Learning curve",
      subtitle = paste0(metric_label, " vs training set size | ",
                        "Gap = overfitting region"),
      x        = "Training samples",
      y        = metric_label
    ) +
    theme_microbiome()

  return(list(plot = p, data = curve_df))
}


# =============================================================================
# SECTION 8 — MODEL STABILITY ANALYSIS
# =============================================================================

#' Assess the stability of feature importance across bootstrap resamples.
#'
#' @param ml_data     Output from prepare_ml_data().
#' @param n_boot      Number of bootstrap iterations. Default = 50.
#' @param top_n       Top features to track. Default = 20.
#' @param n_trees     RF trees per bootstrap. Default = 200.
#' @param seed        Random seed. Default = 42.
#' @return A list: stability data frame and plot.

assess_feature_stability <- function(ml_data,
                                      n_boot  = 50,
                                      top_n   = 20,
                                      n_trees = 200,
                                      seed    = 42) {

  cat("=== Feature importance stability (", n_boot, "bootstraps) ===\n")

  train_df <- ml_data$train
  features <- ml_data$features
  X_train  <- train_df[, features]
  y_train  <- train_df$label

  imp_list <- lapply(seq_len(n_boot), function(i) {
    set.seed(seed + i)
    boot_idx <- sample(nrow(X_train), replace = TRUE)
    X_boot   <- X_train[boot_idx, ]
    y_boot   <- y_train[boot_idx]

    rf <- randomForest(
      x         = X_boot,
      y         = y_boot,
      ntree     = n_trees,
      mtry      = floor(sqrt(ncol(X_boot))),
      importance = TRUE
    )

    imp <- importance(rf, type = 1)  # Mean decrease accuracy
    data.frame(
      feature    = rownames(imp),
      importance = as.numeric(imp[, ncol(imp)]),
      boot       = i,
      stringsAsFactors = FALSE
    )
  })

  all_imp <- bind_rows(imp_list)

  # Stability: coefficient of variation per feature
  stability_df <- all_imp %>%
    group_by(feature) %>%
    summarise(
      mean_imp   = mean(importance),
      sd_imp     = sd(importance),
      cv_imp     = round(sd(importance) / abs(mean(importance)), 4),
      freq_top50 = mean(importance > quantile(importance, 0.5)),
      .groups    = "drop"
    ) %>%
    arrange(desc(mean_imp)) %>%
    slice_head(n = top_n)

  cat("  Top 5 most stable features (low CV):\n")
  print(stability_df %>%
          arrange(cv_imp) %>%
          select(feature, mean_imp, cv_imp) %>%
          head(5))
  cat("\n")

  # Boxplot of importance across bootstraps for top features
  top_feats <- stability_df$feature
  plot_df   <- all_imp %>%
    filter(feature %in% top_feats) %>%
    mutate(feature = factor(feature,
                             levels = rev(stability_df$feature)))

  p <- ggplot(plot_df, aes(x = importance, y = feature)) +
    geom_boxplot(fill = "#3498db", alpha = 0.65, outlier.size = 0.8,
                 outlier.alpha = 0.4, width = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.7) +
    labs(
      title    = paste0("Feature importance stability (", n_boot, " bootstraps)"),
      subtitle = "Wider box = less stable across resamples",
      x        = "Mean decrease in accuracy",
      y        = NULL
    ) +
    theme_microbiome() +
    theme(axis.text.y = element_text(face = "italic", size = 8))

  return(list(
    stability = stability_df,
    plot      = p,
    raw       = all_imp
  ))
}


# =============================================================================
# SECTION 9 — MULTI-CLASS EXTENSION
# =============================================================================

#' Extend binary classification results to multi-class using one-vs-rest.
#'
#' @param ml_data     Output from prepare_ml_data() with n_classes > 2.
#' @param rf_result   Output from train_random_forest().
#' @return A list: per-class ROC curves and macro-averaged AUC.

plot_multiclass_roc <- function(ml_data, rf_result) {

  cat("=== Multi-class ROC analysis (one-vs-rest) ===\n")

  if (ml_data$n_classes <= 2) {
    cat("  Only 2 classes — use plot_roc_curves() instead.\n\n")
    return(NULL)
  }

  y_test  <- ml_data$test$label
  probs   <- rf_result$probabilities

  roc_list <- list()
  auc_vals <- c()

  colours  <- brewer.pal(max(3, ml_data$n_classes), "Set2")[
    seq_len(ml_data$n_classes)
  ]

  for (cls in ml_data$labels) {
    y_bin     <- factor(ifelse(y_test == cls, cls, paste0("not_", cls)),
                         levels = c(paste0("not_", cls), cls))
    prob_cls  <- probs[, cls]
    roc_obj   <- roc(y_bin, prob_cls, quiet = TRUE)
    auc_val   <- round(auc(roc_obj), 4)
    roc_list[[cls]] <- roc_obj
    auc_vals[cls]   <- auc_val
    cat(" ", cls, "(vs rest) AUC:", auc_val, "\n")
  }

  macro_auc <- round(mean(auc_vals), 4)
  cat("  Macro-average AUC:", macro_auc, "\n\n")

  # Plot
  roc_df <- lapply(ml_data$labels, function(cls) {
    r <- roc_list[[cls]]
    data.frame(
      specificity = r$specificities,
      sensitivity = r$sensitivities,
      class       = paste0(cls, " (AUC=", auc_vals[cls], ")"),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity,
                            colour = class)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = "grey70") +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    scale_colour_brewer(palette = "Set2", name = "Class (vs rest)") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    coord_equal() +
    annotate("text", x = 0.75, y = 0.1,
             label = paste0("Macro AUC = ", macro_auc),
             size = 4, fontface = "bold", colour = "grey30") +
    labs(
      title    = "Multi-class ROC — one-vs-rest",
      subtitle = paste0(ml_data$n_classes, " classes | ",
                        "Macro-average AUC = ", macro_auc),
      x        = "1 − Specificity",
      y        = "Sensitivity"
    ) +
    theme_microbiome()

  return(list(
    roc_objects = roc_list,
    auc_values  = auc_vals,
    macro_auc   = macro_auc,
    plot        = p
  ))
}


# =============================================================================
# SECTION 10 — COMPLETE ML WORKFLOW WRAPPER
# =============================================================================

#' Run the complete machine learning classification pipeline.
#'
#' @param ps          A filtered phyloseq object (raw counts).
#' @param group_var   Target variable for classification.
#' @param rank        Taxonomic rank. Default = "Genus".
#' @param test_frac   Test set fraction. Default = 0.2.
#' @param n_trees     RF trees. Default = 500.
#' @param cv_folds    CV folds. Default = 5.
#' @param n_boot      Bootstrap iterations for stability. Default = 50.
#' @param output_dir  Directory for outputs.
#' @return A named list of all results.

run_ml_classification <- function(ps,
                                   group_var  = NULL,
                                   rank       = "Genus",
                                   test_frac  = 0.2,
                                   n_trees    = 500,
                                   cv_folds   = 5,
                                   n_boot     = 50,
                                   output_dir = "ml_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  MACHINE LEARNING CLASSIFICATION PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  results <- list()

  # --- Data preparation -----------------------------------------------------
  ml_data <- prepare_ml_data(
    ps, group_var = group_var, rank = rank,
    test_fraction = test_frac, transform = "clr",
    balance_classes = FALSE
  )
  results$ml_data <- ml_data

  # --- Random Forest --------------------------------------------------------
  cat("--- Training Random Forest ---\n")
  rf_res <- train_random_forest(ml_data, n_trees = n_trees,
                                 cv_folds = cv_folds)
  results$rf <- rf_res

  # --- LASSO ----------------------------------------------------------------
  cat("--- Training LASSO ---\n")
  lasso_res <- train_lasso(ml_data, cv_folds = 10)
  results$lasso <- lasso_res

  # --- Feature importance ---------------------------------------------------
  cat("--- Plot 1: Feature importance ---\n")
  imp_res <- plot_feature_importance(rf_res, lasso_res, top_n = 30, ps = ps)
  results$p_importance <- imp_res$plot
  ggsave(file.path(output_dir, "01_feature_importance.pdf"),
         imp_res$plot, width = 16, height = 12)

  if (!is.null(imp_res$importance$rf)) {
    write.csv(imp_res$importance$rf,
              file.path(output_dir, "rf_feature_importance.csv"),
              row.names = FALSE)
  }
  if (!is.null(imp_res$importance$lasso)) {
    write.csv(imp_res$importance$lasso,
              file.path(output_dir, "lasso_selected_features.csv"),
              row.names = FALSE)
  }

  # --- ROC curves -----------------------------------------------------------
  cat("--- Plot 2: ROC curves ---\n")
  if (ml_data$n_classes == 2) {
    roc_res <- plot_roc_curves(
      model_results = list(RandomForest = rf_res, LASSO = lasso_res),
      ml_data       = ml_data
    )
  } else {
    roc_res <- plot_multiclass_roc(ml_data, rf_res)
  }
  results$p_roc <- roc_res$plot
  ggsave(file.path(output_dir, "02_roc_curves.pdf"),
         roc_res$plot, width = 14, height = 7)

  # --- Confusion matrix -----------------------------------------------------
  cat("--- Plot 3: Confusion matrix ---\n")
  cm_res <- plot_confusion_matrix(rf_res$confusion, model_name = "Random Forest")
  results$p_confusion <- cm_res$plot
  ggsave(file.path(output_dir, "03_confusion_matrix.pdf"),
         cm_res$plot, width = 10, height = 7)

  # --- Learning curve -------------------------------------------------------
  cat("--- Plot 4: Learning curve ---\n")
  lc_res <- plot_learning_curve(ml_data, cv_folds = cv_folds)
  results$p_learning  <- lc_res$plot
  ggsave(file.path(output_dir, "04_learning_curve.pdf"),
         lc_res$plot, width = 9, height = 6)

  # --- Feature stability ----------------------------------------------------
  cat("--- Plot 5: Feature stability ---\n")
  stab_res <- assess_feature_stability(ml_data, n_boot = n_boot)
  results$p_stability <- stab_res$plot
  write.csv(stab_res$stability,
            file.path(output_dir, "feature_stability.csv"),
            row.names = FALSE)
  ggsave(file.path(output_dir, "05_feature_stability.pdf"),
         stab_res$plot, width = 10, height = 8)

  # --- Summary table --------------------------------------------------------
  summary_df <- data.frame(
    Model       = c("Random Forest", "LASSO"),
    Features    = c(length(ml_data$features),
                    lasso_res$n_selected),
    CV_AUC      = c(round(max(rf_res$cv_results[[
      if (ml_data$n_classes == 2) "ROC" else "Accuracy"
    ]]), 4), NA),
    Test_Acc    = c(round(rf_res$confusion$overall["Accuracy"], 4),
                    if (!is.null(lasso_res$confusion))
                      round(lasso_res$confusion$overall["Accuracy"], 4)
                    else NA)
  )

  cat("\n  Summary:\n")
  print(summary_df)
  write.csv(summary_df,
            file.path(output_dir, "model_summary.csv"),
            row.names = FALSE)

  cat("\n", strrep("=", 60), "\n")
  cat("  ML CLASSIFICATION PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  5 PDF files\n")
  cat("  Tables: feature importance, stability, model summary\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_ml_classification(
#   ps         = ps,
#   group_var  = "disease_status",
#   rank       = "Genus",
#   test_frac  = 0.2,
#   n_trees    = 500,
#   cv_folds   = 5,
#   n_boot     = 50,
#   output_dir = "results/ml"
# )

# --- Option B: Individual steps ----------------------------------------------
# ml_data  <- prepare_ml_data(ps, group_var = "disease_status",
#                              rank = "Genus", test_fraction = 0.2)
# rf_res   <- train_random_forest(ml_data, n_trees = 500, cv_folds = 5)
# lasso_res <- train_lasso(ml_data)
# imp_res  <- plot_feature_importance(rf_res, lasso_res, top_n = 30)
# imp_res$plot

# --- Option C: ROC only with external models --------------------------------
# roc_res  <- plot_roc_curves(
#   model_results = list(RF = rf_res, LASSO = lasso_res),
#   ml_data       = ml_data,
#   positive_class = "Disease"
# )
# roc_res$auc_table
# roc_res$plot
