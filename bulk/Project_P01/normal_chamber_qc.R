cell_nature_cols <- list(
  la = "#4DBBD5FF",
  ra = "#E64B35FF",
  ns = "#BEBEBE",
  pair = "#7F7F7F80",
  line = "#7F7F7F",
  heat_low = "#3C5488FF",
  heat_mid = "#F7F7F7",
  heat_high = "#E64B35FF"
)

open_pdf_device <- function(out_file, width, height) {
  pdf(out_file, width = width, height = height, useDingbats = FALSE)
}

safe_row_zscore <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

get_donor_display <- function(sample_id) {
  sub("^[^_]+_", "", as.character(sample_id))
}

draw_pair_segments <- function(dt, x_col, y_col) {
  for (d in unique(dt$donor)) {
    pair_dt <- dt[donor == d]
    if (nrow(pair_dt) != 2) next
    segments(
      x0 = pair_dt[[x_col]][1], y0 = pair_dt[[y_col]][1],
      x1 = pair_dt[[x_col]][2], y1 = pair_dt[[y_col]][2],
      col = cell_nature_cols$pair, lwd = 1.1
    )
  }
}

draw_pair_mid_labels <- function(dt, x_col, y_col, cex = 0.78) {
  xr <- range(dt[[x_col]], na.rm = TRUE)
  yr <- range(dt[[y_col]], na.rm = TRUE)
  y_off <- 0.02 * diff(yr)
  if (!is.finite(y_off) || y_off == 0) y_off <- 0.5

  for (d in unique(dt$donor)) {
    pair_dt <- dt[donor == d]
    if (nrow(pair_dt) != 2) next
    text(
      x = mean(pair_dt[[x_col]], na.rm = TRUE),
      y = mean(pair_dt[[y_col]], na.rm = TRUE) + y_off,
      labels = unique(get_donor_display(pair_dt$sample_id))[1],
      cex = cex,
      font = 2,
      col = cell_nature_cols$line
    )
  }
}

plot_paired_scatter <- function(dt, x_col, y_col, xlab, ylab, main, out_file, width = 9.5, height = 8.0, label_mode = c("sample", "donor_mid", "none")) {
  label_mode <- match.arg(label_mode)
  open_pdf_device(out_file, width = width, height = height)
  cols <- ifelse(dt$chamber == "LA", cell_nature_cols$la, cell_nature_cols$ra)
  plot(dt[[x_col]], dt[[y_col]], pch = 19, col = cols, xlab = xlab, ylab = ylab, main = main)
  draw_pair_segments(dt, x_col, y_col)
  if (label_mode == "sample") {
    text(dt[[x_col]], dt[[y_col]], labels = dt$sample_id, pos = 3, cex = 0.75)
  } else if (label_mode == "donor_mid") {
    draw_pair_mid_labels(dt, x_col, y_col)
  }
  legend("topright", legend = c("LA", "RA"), col = c(cell_nature_cols$la, cell_nature_cols$ra), pch = 19, bty = "n")
  dev.off()
}

plot_figure1a <- function(pca_dt, out_file) {
  plot_paired_scatter(
    dt = pca_dt,
    x_col = "PC1",
    y_col = "PC2",
    xlab = "PC1",
    ylab = "PC2",
    main = "Figure 1A: GSE112339 PCA (LA vs RA)",
    out_file = out_file,
    label_mode = "donor_mid"
  )
}

plot_correlation_heatmap <- function(corr_mat, out_file) {
  open_pdf_device(out_file, width = 9.5, height = 8.5)
  ord <- hclust(as.dist(1 - corr_mat))$order
  pal <- colorRampPalette(c("white", cell_nature_cols$heat_low, cell_nature_cols$heat_high))(100)
  heatmap(corr_mat[ord, ord], Rowv = NA, Colv = NA, scale = "none", margins = c(8, 8), col = pal)
  dev.off()
}

write_pca_chamber_lm <- function(pca_dt, out_file_txt, out_file_tsv) {
  df_pca <- as.data.frame(pca_dt)
  df_pca$chamber <- factor(df_pca$chamber, levels = c("RA", "LA"))
  df_pca$donor <- factor(df_pca$donor)

  fit_pc1 <- lm(PC1 ~ chamber + donor, data = df_pca)
  fit_pc2 <- lm(PC2 ~ chamber + donor, data = df_pca)

  txt <- c(
    "summary(lm(PC1 ~ chamber + donor, data = df_pca))",
    capture.output(summary(fit_pc1)),
    "",
    "summary(lm(PC2 ~ chamber + donor, data = df_pca))",
    capture.output(summary(fit_pc2))
  )
  writeLines(txt, con = out_file_txt)

  tidy_one <- function(fit, response_name) {
    cf <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    setnames(cf, c("term", "estimate", "std_error", "t_value", "p_value"))
    cf[, response := response_name]
    cf
  }

  lm_dt <- rbindlist(list(tidy_one(fit_pc1, "PC1"), tidy_one(fit_pc2, "PC2")), fill = TRUE)
  fwrite(lm_dt[, .(response, term, estimate, std_error, t_value, p_value)], out_file_tsv, sep = "\t")
}

plot_logcpm_mean_variance <- function(log_cpm_mat, out_file) {
  gene_mean <- rowMeans(log_cpm_mat)
  gene_var <- apply(log_cpm_mat, 1, var)
  ord <- order(gene_mean)
  smooth_xy <- lowess(gene_mean[ord], gene_var[ord], f = 0.2)

  open_pdf_device(out_file, width = 9.5, height = 8.0)
  plot(
    gene_mean, gene_var,
    pch = 16, cex = 0.35,
    col = adjustcolor(cell_nature_cols$line, alpha.f = 0.28),
    xlab = "rowMeans(logCPM)",
    ylab = "apply(logCPM, 1, var)",
    main = "Gene Mean-Variance Diagnostic"
  )
  lines(smooth_xy$x, smooth_xy$y, col = cell_nature_cols$ra, lwd = 2.2)
  dev.off()
}

plot_figure1b <- function(gene_out, expr_mat, meta, out_file, gene_fdr_cutoff = 0.05, top_each = 20, anchor_logfc = 0.5) {
  volc <- copy(gene_out)
  volc[, neglog10fdr := -log10(pmax(fdr, 1e-300))]
  volc[, sig_class := fifelse(fdr < gene_fdr_cutoff & logFC > 0, "LA-high", fifelse(fdr < gene_fdr_cutoff & logFC < 0, "RA-high", "NS"))]
  label_each <- 6L

  la_top <- volc[fdr < gene_fdr_cutoff & logFC > anchor_logfc][order(-moderated_t)]
  ra_top <- volc[fdr < gene_fdr_cutoff & logFC < -anchor_logfc][order(moderated_t)]
  if (nrow(la_top) < top_each) la_top <- volc[logFC > 0][order(-moderated_t)]
  if (nrow(ra_top) < top_each) ra_top <- volc[logFC < 0][order(moderated_t)]
  la_top <- head(la_top, top_each)
  ra_top <- head(ra_top, top_each)
  sel_keys <- unique(c(la_top$gene_key, ra_top$gene_key))

  hm <- matrix(numeric(0), nrow = 0, ncol = 0)
  row_labels <- character(0)
  meta_ord <- copy(meta)
  if (length(sel_keys) > 0) {
    meta_ord <- copy(meta)
    meta_ord[, chamber_plot := factor(chamber, levels = c("RA", "LA"))]
    setorder(meta_ord, donor, chamber_plot, sample_id)
    hm <- expr_mat[sel_keys, meta_ord$sample_id, drop = FALSE]
    hm <- safe_row_zscore(hm)
    row_dt <- volc[match(rownames(hm), gene_key)]
    row_labels <- fifelse(!is.na(row_dt$gene_symbol) & nzchar(row_dt$gene_symbol), row_dt$gene_symbol, row_dt$ensembl_id)
    row_labels <- make.unique(as.character(row_labels))
    row_ord <- order(row_dt$logFC, decreasing = TRUE)
    hm <- hm[row_ord, , drop = FALSE]
    row_labels <- row_labels[row_ord]
  }

  open_pdf_device(out_file, width = 16, height = 8)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 1.25))

  par(mar = c(5, 5, 3, 1))
  pt_col <- ifelse(volc$sig_class == "LA-high", cell_nature_cols$la, ifelse(volc$sig_class == "RA-high", cell_nature_cols$ra, cell_nature_cols$ns))
  plot(volc$logFC, volc$neglog10fdr, pch = 16, cex = 0.55, col = pt_col, xlab = "logFC (LA relative to RA)", ylab = "-log10(FDR)", main = "Figure 1B (left): Volcano")
  abline(v = c(-0.5, 0.5), col = "#636363", lty = 2)
  abline(h = -log10(gene_fdr_cutoff), col = "#636363", lty = 2)
  la_label <- volc[fdr < gene_fdr_cutoff & logFC > anchor_logfc][order(-moderated_t)]
  if (nrow(la_label) < label_each) la_label <- volc[logFC > 0][order(-moderated_t)]
  la_label <- head(la_label, label_each)
  ra_label <- volc[fdr < gene_fdr_cutoff & logFC < -anchor_logfc][order(moderated_t)]
  if (nrow(ra_label) < label_each) ra_label <- volc[logFC < 0][order(moderated_t)]
  ra_label <- head(ra_label, label_each)
  top_label <- unique(rbindlist(list(la_label, ra_label), fill = TRUE), by = "gene_key")
  if (nrow(top_label) > 0) {
    lbl <- fifelse(!is.na(top_label$gene_symbol) & nzchar(top_label$gene_symbol), top_label$gene_symbol, top_label$ensembl_id)
    text(top_label$logFC, top_label$neglog10fdr, labels = lbl, pos = 3, cex = 0.65)
  }
  legend("topright", legend = c("LA-high (FDR<0.05)", "RA-high (FDR<0.05)", "NS"), col = c(cell_nature_cols$la, cell_nature_cols$ra, cell_nature_cols$ns), pch = 16, bty = "n", cex = 0.85)

  par(mar = c(10, 9, 4.5, 2))
  if (nrow(hm) == 0 || ncol(hm) == 0) {
    plot.new()
    title(main = "Figure 1B (right): Top chamber genes heatmap")
    text(0.5, 0.5, "No genes passed the current selection rule")
  } else {
    nc <- ncol(hm)
    nr <- nrow(hm)
    pal <- colorRampPalette(c(cell_nature_cols$heat_low, cell_nature_cols$heat_mid, cell_nature_cols$heat_high))(100)
    lim <- max(abs(hm))
    image(x = seq_len(nc), y = seq_len(nr), z = t(hm[nr:1, , drop = FALSE]), col = pal, zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "", main = "Figure 1B (right): Top chamber genes heatmap")
    axis(1, at = seq_len(nc), labels = meta_ord$sample_id, las = 2, cex.axis = if (nc > 20) 0.45 else 0.7)
    axis(2, at = seq_len(nr), labels = rev(row_labels), las = 2, cex.axis = if (nr > 35) 0.4 else 0.6)
    box()
    donor_breaks <- which(meta_ord$donor[-1] != meta_ord$donor[-nrow(meta_ord)])
    if (length(donor_breaks) > 0) {
      abline(v = donor_breaks + 0.5, lwd = 1.4, col = "#4D4D4D")
    }
    donor_blocks <- split(seq_len(nc), meta_ord$donor)
    donor_labels <- vapply(donor_blocks, function(idx) get_donor_display(meta_ord$sample_id[idx[1]]), character(1))
    usr <- par("usr")
    y_top <- usr[4]
    y_tick <- y_top + 0.03 * (usr[4] - usr[3])
    y_lab <- y_top + 0.08 * (usr[4] - usr[3])
    for (i in seq_along(donor_blocks)) {
      idx <- donor_blocks[[i]]
      x0 <- min(idx) - 0.35
      x1 <- max(idx) + 0.35
      segments(x0 = x0, y0 = y_top, x1 = x0, y1 = y_tick, xpd = NA, lwd = 1.2, col = cell_nature_cols$line)
      segments(x0 = x1, y0 = y_top, x1 = x1, y1 = y_tick, xpd = NA, lwd = 1.2, col = cell_nature_cols$line)
      segments(x0 = x0, y0 = y_tick, x1 = x1, y1 = y_tick, xpd = NA, lwd = 1.2, col = cell_nature_cols$line)
      text(x = mean(idx), y = y_lab, labels = donor_labels[i], xpd = NA, cex = 0.78, font = 2, col = cell_nature_cols$line)
    }
  }
  dev.off()
}

plot_figure1b_compact <- function(gene_out, expr_mat, meta, out_file, gene_fdr_cutoff = 0.05, top_each = 6, anchor_logfc = 0.5) {
  volc <- copy(gene_out)
  volc[, neglog10fdr := -log10(pmax(fdr, 1e-300))]
  volc[, sig_class := fifelse(fdr < gene_fdr_cutoff & logFC > 0, "LA-high", fifelse(fdr < gene_fdr_cutoff & logFC < 0, "RA-high", "NS"))]
  label_each <- 4L

  la_top <- volc[fdr < gene_fdr_cutoff & logFC > anchor_logfc][order(-moderated_t)]
  ra_top <- volc[fdr < gene_fdr_cutoff & logFC < -anchor_logfc][order(moderated_t)]
  if (nrow(la_top) < top_each) la_top <- volc[logFC > 0][order(-moderated_t)]
  if (nrow(ra_top) < top_each) ra_top <- volc[logFC < 0][order(moderated_t)]
  la_top <- head(la_top, top_each)
  ra_top <- head(ra_top, top_each)
  sel_keys <- unique(c(la_top$gene_key, ra_top$gene_key))

  hm <- matrix(numeric(0), nrow = 0, ncol = 0)
  row_labels <- character(0)
  meta_ord <- copy(meta)
  if (length(sel_keys) > 0) {
    meta_ord <- copy(meta)
    meta_ord[, chamber_plot := factor(chamber, levels = c("RA", "LA"))]
    setorder(meta_ord, donor, chamber_plot, sample_id)
    hm <- expr_mat[sel_keys, meta_ord$sample_id, drop = FALSE]
    hm <- safe_row_zscore(hm)
    row_dt <- volc[match(rownames(hm), gene_key)]
    row_labels <- fifelse(!is.na(row_dt$gene_symbol) & nzchar(row_dt$gene_symbol), row_dt$gene_symbol, row_dt$ensembl_id)
    row_labels <- make.unique(as.character(row_labels))
    row_ord <- order(row_dt$logFC, decreasing = TRUE)
    hm <- hm[row_ord, , drop = FALSE]
    row_labels <- row_labels[row_ord]
  }

  open_pdf_device(out_file, width = 12, height = 5.8)
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.82, 1.38))

  par(mar = c(4.2, 4.6, 2.6, 0.8))
  pt_col <- ifelse(volc$sig_class == "LA-high", cell_nature_cols$la, ifelse(volc$sig_class == "RA-high", cell_nature_cols$ra, cell_nature_cols$ns))
  plot(
    volc$logFC, volc$neglog10fdr,
    pch = 16, cex = 0.65, col = pt_col,
    xlab = "logFC (LA vs RA)", ylab = "-log10(FDR)",
    main = "Volcano"
  )
  abline(v = c(-0.5, 0.5), col = "#636363", lty = 2)
  abline(h = -log10(gene_fdr_cutoff), col = "#636363", lty = 2)
  la_label <- volc[fdr < gene_fdr_cutoff & logFC > anchor_logfc][order(-moderated_t)]
  if (nrow(la_label) < label_each) la_label <- volc[logFC > 0][order(-moderated_t)]
  la_label <- head(la_label, label_each)
  ra_label <- volc[fdr < gene_fdr_cutoff & logFC < -anchor_logfc][order(moderated_t)]
  if (nrow(ra_label) < label_each) ra_label <- volc[logFC < 0][order(moderated_t)]
  ra_label <- head(ra_label, label_each)
  top_label <- unique(rbindlist(list(la_label, ra_label), fill = TRUE), by = "gene_key")
  if (nrow(top_label) > 0) {
    lbl <- fifelse(!is.na(top_label$gene_symbol) & nzchar(top_label$gene_symbol), top_label$gene_symbol, top_label$ensembl_id)
    text(top_label$logFC, top_label$neglog10fdr, labels = lbl, pos = 3, cex = 0.72)
  }
  legend(
    "topright",
    legend = c("LA-high", "RA-high", "NS"),
    col = c(cell_nature_cols$la, cell_nature_cols$ra, cell_nature_cols$ns),
    pch = 16, bty = "n", cex = 0.82
  )

  par(mar = c(2.5, 7.8, 5.2, 2.4))
  if (nrow(hm) == 0 || ncol(hm) == 0) {
    plot.new()
    title(main = "Top paired chamber genes")
    text(0.5, 0.5, "No genes passed the current selection rule")
  } else {
    nc <- ncol(hm)
    nr <- nrow(hm)
    pal <- colorRampPalette(c(cell_nature_cols$heat_low, cell_nature_cols$heat_mid, cell_nature_cols$heat_high))(100)
    lim <- max(abs(hm))
    image(
      x = seq_len(nc), y = seq_len(nr), z = t(hm[nr:1, , drop = FALSE]),
      col = pal, zlim = c(-lim, lim), axes = FALSE,
      xlab = "", ylab = "", main = "Top paired chamber genes"
    )
    axis(2, at = seq_len(nr), labels = rev(row_labels), las = 2, cex.axis = 0.8)
    box()

    donor_breaks <- which(meta_ord$donor[-1] != meta_ord$donor[-nrow(meta_ord)])
    if (length(donor_breaks) > 0) {
      abline(v = donor_breaks + 0.5, lwd = 1.1, col = "#4D4D4D")
    }

    usr <- par("usr")
    y_top <- usr[4]
    y_bar0 <- y_top + 0.10 * (usr[4] - usr[3])
    y_bar1 <- y_top + 0.22 * (usr[4] - usr[3])
    y_bracket <- y_top + 0.34 * (usr[4] - usr[3])
    y_donor <- y_top + 0.46 * (usr[4] - usr[3])

    chamber_cols <- c(RA = cell_nature_cols$ra, LA = cell_nature_cols$la)
    for (i in seq_len(nc)) {
      rect(i - 0.5, y_bar0, i + 0.5, y_bar1, col = chamber_cols[as.character(meta_ord$chamber[i])], border = NA, xpd = NA)
    }
    text(x = 0.15, y = (y_bar0 + y_bar1) / 2, labels = "Chamber", adj = c(0, 0.5), xpd = NA, cex = 0.72, col = cell_nature_cols$line)

    donor_blocks <- split(seq_len(nc), meta_ord$donor)
    donor_labels <- vapply(donor_blocks, function(idx) get_donor_display(meta_ord$sample_id[idx[1]]), character(1))
    for (i in seq_along(donor_blocks)) {
      idx <- donor_blocks[[i]]
      x0 <- min(idx) - 0.35
      x1 <- max(idx) + 0.35
      segments(x0 = x0, y0 = y_top, x1 = x0, y1 = y_bracket, xpd = NA, lwd = 1.1, col = cell_nature_cols$line)
      segments(x0 = x1, y0 = y_top, x1 = x1, y1 = y_bracket, xpd = NA, lwd = 1.1, col = cell_nature_cols$line)
      segments(x0 = x0, y0 = y_bracket, x1 = x1, y1 = y_bracket, xpd = NA, lwd = 1.1, col = cell_nature_cols$line)
      text(x = mean(idx), y = y_donor, labels = donor_labels[i], xpd = NA, cex = 0.74, font = 2, col = cell_nature_cols$line)
    }

    legend(
      "bottom",
      inset = c(0, -0.18),
      legend = c("RA", "LA"),
      fill = c(cell_nature_cols$ra, cell_nature_cols$la),
      border = NA,
      bty = "n",
      horiz = TRUE,
      xpd = NA,
      cex = 0.8
    )
  }
  dev.off()
}

run_qc_outputs <- function(dge, meta, qc_dir, fig_dir) {
  qc_lib <- data.table(
    sample_id = meta$sample_id,
    donor = as.character(meta$donor),
    chamber = as.character(meta$chamber),
    library_size = colSums(dge$counts),
    detected_genes = colSums(dge$counts > 0)
  )
  fwrite(qc_lib, file.path(qc_dir, "qc_library_and_detected_genes.tsv"), sep = "\t")

  log_cpm_all <- cpm(dge, log = TRUE, prior.count = 1)
  corr <- cor(log_cpm_all, method = "spearman")
  corr_dt <- as.data.table(as.table(corr))
  setnames(corr_dt, c("sample_1", "sample_2", "spearman_corr"))
  fwrite(corr_dt, file.path(qc_dir, "qc_sample_correlation.tsv"), sep = "\t")
  plot_correlation_heatmap(corr, file.path(qc_dir, "qc_sample_correlation_heatmap.pdf"))

  pca <- prcomp(t(log_cpm_all), scale. = TRUE)
  pca_dt <- data.table(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  pca_dt <- merge(pca_dt, meta[, .(sample_id, donor, chamber)], by = "sample_id", all.x = TRUE)
  fwrite(pca_dt, file.path(qc_dir, "qc_pca.tsv"), sep = "\t")
  write_pca_chamber_lm(pca_dt, file.path(qc_dir, "qc_pca_chamber_lm.txt"), file.path(qc_dir, "qc_pca_chamber_lm_coefficients.tsv"))
  plot_figure1a(pca_dt, file.path(fig_dir, "figure1A_gse112339_pca_la_vs_ra.pdf"))
  plot_paired_scatter(pca_dt, "PC1", "PC2", "PC1", "PC2", "PCA (logCPM)", file.path(qc_dir, "qc_pca.pdf"))
  plot_logcpm_mean_variance(log_cpm_all, file.path(qc_dir, "qc_logcpm_mean_variance.pdf"))

  mds <- plotMDS(log_cpm_all, plot = FALSE)
  mds_dt <- data.table(sample_id = colnames(log_cpm_all), MDS1 = mds$x, MDS2 = mds$y)
  mds_dt <- merge(mds_dt, meta[, .(sample_id, donor, chamber)], by = "sample_id", all.x = TRUE)
  fwrite(mds_dt, file.path(qc_dir, "qc_mds.tsv"), sep = "\t")
  plot_paired_scatter(mds_dt, "MDS1", "MDS2", "MDS1", "MDS2", "MDS", file.path(qc_dir, "qc_mds.pdf"))

  list(log_cpm_all = log_cpm_all)
}
