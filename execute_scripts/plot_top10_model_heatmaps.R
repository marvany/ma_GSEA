#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) >= 1) args[[1]] else "rachel_GSEA/Results_260212"
out_dir <- if (length(args) >= 2) args[[2]] else file.path(results_root, "All.Plots", "top10_model_single_column_heatmaps")
top_n <- if (length(args) >= 3) as.integer(args[[3]]) else 10L

if (is.na(top_n) || top_n <= 0) stop("top_n must be a positive integer")
if (!dir.exists(results_root)) stop("Results root does not exist: ", results_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

list_result_files <- function(root) {
  files <- list.files(
    root,
    pattern = "\\.genesets\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )
  files[grepl("result_textFiles_individual", files, fixed = TRUE)]
}

extract_meta <- function(file, root) {
  rel <- sub(paste0("^", normalizePath(root, winslash = "/", mustWork = TRUE), "/"), "", normalizePath(file, winslash = "/", mustWork = TRUE))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]

  # Expected structure:
  # pvalue_xx/GWAS/model_ID/geneset/result_textFiles_individual/file.tsv
  if (length(parts) < 6L) return(NULL)
  list(
    pvalue = parts[[1]],
    gwas = parts[[2]],
    model_ID = parts[[3]],
    geneset = parts[[4]]
  )
}

clean_pathway <- function(x) {
  out <- gsub("[._]+", " ", x)
  out <- trimws(out)
  tools::toTitleCase(tolower(out))
}

wrap_text <- function(x, width = 34L) {
  vapply(
    x,
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    FUN.VALUE = character(1)
  )
}

read_one_result <- function(file, root) {
  meta <- extract_meta(file, root)
  if (is.null(meta)) return(NULL)

  dt <- tryCatch(fread(file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)

  if (!("Reference" %in% names(dt)) || !("pval" %in% names(dt))) return(NULL)

  dt <- dt[, .(Reference, pval)]
  dt <- dt[is.finite(pval) & pval > 0]
  if (nrow(dt) == 0L) return(NULL)

  dt[, `:=`(
    pvalue = meta$pvalue,
    gwas = meta$gwas,
    model_ID = meta$model_ID,
    geneset = meta$geneset,
    pathway = clean_pathway(Reference),
    neg_log10p = -log10(pval)
  )]

  dt[, .(pvalue, gwas, model_ID, geneset, pathway, pval, neg_log10p)]
}

build_top_table <- function(root, top_n = 10L) {
  files <- list_result_files(root)
  if (length(files) == 0L) stop("No result files found under: ", root)

  dt_list <- lapply(files, read_one_result, root = root)
  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  if (nrow(dt) == 0L) stop("Result files were found, but no usable rows were parsed.")

  setorder(dt, pvalue, gwas, model_ID, pval)
  dt[, head(.SD, top_n), by = .(pvalue, gwas, model_ID)]
}

plot_one_combo <- function(dt_combo, combo_pvalue, combo_gwas, out_dir) {
  dtp <- copy(dt_combo)
  if (nrow(dtp) == 0L) return(invisible(NULL))
  n_models <- uniqueN(dtp$model_ID)
  n_cols <- min(3L, n_models)
  n_rows <- ceiling(n_models / n_cols)

  setorder(dtp, model_ID, pval)

  dtp[, tile_x := " "]
  dtp[, pathway_display := wrap_text(pathway, width = 34L)]
  # Unique y labels per panel with stable within-model ordering by significance.
  dtp[, pathway_panel := paste(model_ID, pathway_display, sep = "|||")]
  dtp[, pathway_panel := factor(pathway_panel, levels = rev(unique(pathway_panel)))]

  max_fill <- max(dtp$neg_log10p, na.rm = TRUE)
  fill_limits <- c(0, ceiling(max_fill))
  if (fill_limits[2] <= 0) fill_limits[2] <- max_fill

  p <- ggplot(dtp, aes(x = tile_x, y = pathway_panel, fill = neg_log10p)) +
    geom_tile(color = "black", linewidth = 0.3, width = 0.5, height = 0.92) +
    facet_wrap(~ model_ID, ncol = n_cols, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_fill_gradientn(
      colors = c("#F5F5F5", "#F0E442", "#E69F00"),
      name = expression(-log[10](italic(P))),
      na.value = "gray50",
      limits = fill_limits,
      breaks = pretty_breaks(n = 5)
    ) +
    labs(
      title = sprintf("%s | %s | Top %d pathways per model", combo_gwas, combo_pvalue, top_n),
      x = NULL,
      y = "Pathway"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.6, "lines"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  file_stub <- sprintf("top10_pathway_heatmap_%s_%s", combo_pvalue, combo_gwas)
  width_in <- max(16, n_cols * 6)
  height_in <- max(10, n_rows * 4.2)

  ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), p, width = width_in, height = height_in, units = "in")
  ggsave(file.path(out_dir, paste0(file_stub, ".png")), p, width = width_in, height = height_in, units = "in", dpi = 300)

  invisible(file_stub)
}

main <- function() {
  dt_top <- build_top_table(results_root, top_n = top_n)
  combos <- unique(dt_top[, .(pvalue, gwas)])
  setorder(combos, pvalue, gwas)

  produced <- character()
  for (i in seq_len(nrow(combos))) {
    pv <- combos$pvalue[[i]]
    gw <- combos$gwas[[i]]
    dt_combo <- dt_top[pvalue == pv & gwas == gw]
    out <- plot_one_combo(dt_combo, pv, gw, out_dir)
    if (!is.null(out)) produced <- c(produced, out)
  }

  summary_file <- file.path(out_dir, "top10_heatmap_summary.tsv")
  fwrite(dt_top, summary_file, sep = "\t")

  message("Done. Generated ", length(produced), " plot(s) in: ", out_dir)
  message("Summary table: ", summary_file)
}

main()
