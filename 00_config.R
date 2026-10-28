# =============================================================================
# 00_config.R  -  Configuracion y funciones comunes (TFM neuroblastoma)
# Se carga desde los demas scripts con source(). NO hace falta ejecutarlo solo.
# =============================================================================
# ROOT ya viene definido por cada script (variable de entorno TFM_ROOT o ruta por defecto)
if (!exists("ROOT")) ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
DATA <- file.path(ROOT, "Datos")          # Datos/18271, Datos/212877, Datos/62564 (y opcional Datos/49711)
RES  <- file.path(ROOT, "resultados")     # todo lo que generan los scripts va aqui
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

options(stringsAsFactors = FALSE, width = 140, warn = 1, scipen = 999)
SEED <- 1234
set.seed(SEED)

CANDIDATES <- c("ARPP21", "DLG2", "AXL", "MYCN")   # genes candidatos del TFM
FDR_CUT <- 0.05                                     # umbral de FDR (BH)
LFC_CUT <- 1                                        # umbral |log2FC| para listas de DEGs

# ---- utilidades ------------------------------------------------------------
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Faltan paquetes: ", paste(miss, collapse = ", "),
                         "\nEjecuta primero 00_instalar_paquetes.R", call. = FALSE)
  invisible(TRUE)
}
has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

.LOG <- new.env()
start_log <- function(name) {
  while (sink.number() > 0) sink()
  .LOG$name <- name
  .LOG$dir  <- file.path(RES, name); dir.create(.LOG$dir, showWarnings = FALSE, recursive = TRUE)
  .LOG$file <- file.path(.LOG$dir, paste0("LOG_", name, ".txt"))
  .LOG$sum  <- file.path(.LOG$dir, paste0("RESUMEN_", name, ".txt"))
  if (file.exists(.LOG$sum)) file.remove(.LOG$sum)
  sink(.LOG$file, split = TRUE)
  cat("== ", name, " == ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat(R.version.string, "\n")
  invisible(.LOG$dir)
}
OUT <- function(...) file.path(.LOG$dir, ...)        # ruta de salida del script en curso
# S(): imprime y ademas guarda la linea en RESUMEN_<name>.txt (lo que me tienes que pegar)
S <- function(...) { txt <- paste0(..., collapse = ""); cat(txt, "\n", sep = ""); cat(txt, "\n", file = .LOG$sum, append = TRUE, sep = "") }
SH <- function(title) S("\n### ", title)
STAB <- function(df, digits = 4) {                    # tabla -> texto en el resumen
  txt <- capture.output(print(format(as.data.frame(df), digits = digits, justify = "left"), row.names = TRUE))
  for (l in txt) S(l)
}
end_log <- function() {
  si <- capture.output(sessionInfo()); writeLines(si, OUT("sessionInfo.txt"))
  cat("\nTerminado. Resultados en: ", .LOG$dir, "\n", sep = "")
  cat("Pega aqui el contenido de: ", .LOG$sum, "\n", sep = "")
  while (sink.number() > 0) sink()
}

save_png <- function(file, expr, width = 1800, height = 1300, res = 200) {
  png(file, width = width, height = height, res = res)
  on.exit(dev.off())
  r <- tryCatch(expr, error = function(e) { message("  (grafico omitido: ", conditionMessage(e), ")"); NULL })
  if (inherits(r, "ggplot")) print(r)
  invisible(NULL)
}

# ---- lector generico de series matrix (GEO) --------------------------------
# Devuelve data.frame con una fila por muestra: title, geo_accession, source y una columna
# por cada clave "clave: valor" de las lineas !Sample_characteristics_ch1 (aunque esten desordenadas).
parse_series_matrix <- function(file) {
  L <- readLines(file, warn = FALSE)
  L <- L[startsWith(L, "!Sample_")]
  row <- function(key) lapply(L[startsWith(L, key)], function(x) gsub('^"|"$', "", strsplit(x, "\t")[[1]][-1]))
  one <- function(key) { r <- row(key); if (length(r)) r[[1]] else NA }
  title <- one("!Sample_title"); n <- length(title)
  out <- data.frame(title = title, geo_accession = one("!Sample_geo_accession"),
                    source = one("!Sample_source_name_ch1"), stringsAsFactors = FALSE)
  ch <- row("!Sample_characteristics_ch1")
  if (length(ch)) {
    ch <- do.call(cbind, ch)                       # n muestras x k filas
    kv <- apply(ch, 1, function(r) {
      r <- r[nzchar(r) & grepl(":", r)]
      k <- trimws(sub(":.*$", "", r)); v <- trimws(sub("^[^:]*:", "", r)); setNames(v, k)
    }, simplify = FALSE)
    keys <- unique(unlist(lapply(kv, names)))
    for (k in keys) out[[make.names(k)]] <- vapply(kv, function(x) if (k %in% names(x)) x[[k]] else NA_character_, "")
  }
  out
}

# ---- tablas de resultados ---------------------------------------------------
deg_counts <- function(tt, fdr = FDR_CUT, lfc = LFC_CUT) {
  c(testados = nrow(tt),
    FDR = sum(tt$adj.P.Val < fdr),
    FDR_y_LFC = sum(tt$adj.P.Val < fdr & abs(tt$logFC) > lfc),
    up = sum(tt$adj.P.Val < fdr & tt$logFC > lfc),
    down = sum(tt$adj.P.Val < fdr & tt$logFC < -lfc))
}
cand_table <- function(tt, genes = CANDIDATES) {
  g <- genes[genes %in% rownames(tt)]
  if (!length(g)) return(NULL)
  d <- tt[g, intersect(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), colnames(tt)), drop = FALSE]
  d[setdiff(genes, g), ] <- NA
  d[genes, , drop = FALSE]
}

volcano_plot <- function(tt, title, file, fdr = FDR_CUT, lfc = LFC_CUT, label = CANDIDATES, xlab = "log2 fold change") {
  need("ggplot2")
  d <- tt; d$gene <- rownames(d)
  d$cat <- "NS"; d$cat[d$adj.P.Val < fdr & d$logFC >  lfc] <- "Up"; d$cat[d$adj.P.Val < fdr & d$logFC < -lfc] <- "Down"
  d$ml <- -log10(pmax(d$adj.P.Val, 1e-300))
  lab <- d[d$gene %in% label, ]
  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, ml, colour = cat)) +
    ggplot2::geom_point(alpha = .5, size = 1) +
    ggplot2::scale_colour_manual(values = c(Down = "#2C7BB6", NS = "grey70", Up = "#D7191C")) +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed") +
    ggplot2::labs(title = title, x = xlab, y = "-log10(FDR)", colour = "") + ggplot2::theme_minimal(base_size = 13)
  if (nrow(lab) && has("ggrepel")) p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), colour = "black", max.overlaps = Inf)
  save_png(file, p, 1800, 1300, 200)
}

pca_plot <- function(mat, group, file, title, shape = NULL, legend_title = "Grupo", subtitle = NULL) {
  need("ggplot2")
  pc <- prcomp(t(mat), scale. = TRUE); ve <- summary(pc)$importance[2, ] * 100
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = group, sample = rownames(pc$x))
  cm <- function(x) sub("-", "\u2212", format(x, decimal.mark = ",", trim = TRUE))
  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = group)) + ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_x_continuous(labels = cm) + ggplot2::scale_y_continuous(labels = cm) +
    ggplot2::labs(title = title, subtitle = subtitle, colour = legend_title,
                  x = paste0("PC1 (", cm(round(ve[1], 1)), " %)"), y = paste0("PC2 (", cm(round(ve[2], 1)), " %)")) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13), plot.title.position = "plot",
                   plot.subtitle = ggplot2::element_text(colour = "grey30", size = 10.5), panel.grid.minor = ggplot2::element_blank())
  save_png(file, p, 1600, 1200, 200)
  invisible(ve)
}

heatmap_png <- function(mat, annot, file, title, scale = "row", ...) {
  if (!has("pheatmap")) { message("  (pheatmap no instalado: heatmap omitido)"); return(invisible()) }
  png(file, 1700, 1500, res = 200); on.exit(dev.off())
  tryCatch(pheatmap::pheatmap(mat, scale = scale, annotation_col = annot, main = title, show_rownames = FALSE, ...),
           error = function(e) message("  (heatmap omitido: ", conditionMessage(e), ")"))
}

# ---- enriquecimiento funcional (ORA + GSEA) -- comun a los tres datasets ----
# tt: data.frame de topTable con rownames = SIMBOLOS unicos y columnas logFC, t, adj.P.Val.
# ORA: genes up / down por separado, universo = genes realmente testados. GSEA: ranking por t moderado.
run_enrichment <- function(tt, prefix, outdir = .LOG$dir, fdr = FDR_CUT, lfc = LFC_CUT, kegg = TRUE, gsea = TRUE) {
  if (!all(vapply(c("clusterProfiler", "org.Hs.eg.db"), has, logical(1)))) { S("  (enriquecimiento omitido: falta clusterProfiler/org.Hs.eg.db)"); return(invisible(NULL)) }
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  set.seed(SEED)
  uni <- rownames(tt)
  res <- list()
  plot_top <- function(obj, file, title, n = 15) {
    if (!is.null(obj) && nrow(as.data.frame(obj)) > 0 && has("enrichplot"))
      save_png(file, suppressMessages(enrichplot::dotplot(obj, showCategory = n) + ggplot2::ggtitle(title)), 2000, 1600, 200)
  }
  for (dir in c("up", "down")) {
    g <- rownames(tt)[tt$adj.P.Val < fdr & (if (dir == "up") tt$logFC > lfc else tt$logFC < -lfc)]
    S(sprintf("  [%s] ORA %s: %d genes", prefix, dir, length(g)))
    if (length(g) < 10) next
    ego <- tryCatch(enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = uni,
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE),
                    error = function(e) { S("   enrichGO fallo: ", conditionMessage(e)); NULL })
    if (!is.null(ego)) {
      d <- as.data.frame(ego); write.csv(d, file.path(outdir, sprintf("%s_ORA_GOBP_%s.csv", prefix, dir)), row.names = FALSE)
      ego2 <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) ego)
      write.csv(as.data.frame(ego2), file.path(outdir, sprintf("%s_ORA_GOBP_%s_simplificado.csv", prefix, dir)), row.names = FALSE)
      S(sprintf("   GO BP significativos (q<0.05): %d  (tras simplify: %d)", nrow(d), nrow(as.data.frame(ego2))))
      if (nrow(d)) { top <- head(as.data.frame(ego2)[, c("Description", "p.adjust")], 8); for (i in seq_len(nrow(top))) S(sprintf("     %-70s q=%.2g", substr(top$Description[i], 1, 70), top$p.adjust[i])) }
      plot_top(ego2, file.path(outdir, sprintf("%s_ORA_GOBP_%s.png", prefix, dir)), sprintf("%s: GO BP (%s)", prefix, dir))
    }
    if (kegg) {
      ent <- suppressWarnings(tryCatch(bitr(g, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db), error = function(e) NULL))
      entu <- suppressWarnings(tryCatch(bitr(uni, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db), error = function(e) NULL))
      if (!is.null(ent) && !is.null(entu)) {
        ek <- tryCatch(enrichKEGG(unique(ent$ENTREZID), organism = "hsa", universe = unique(entu$ENTREZID), pAdjustMethod = "BH", pvalueCutoff = 0.05),
                       error = function(e) { S("   enrichKEGG no disponible (sin internet/KEGG): ", substr(conditionMessage(e), 1, 80)); NULL })
        if (!is.null(ek)) { write.csv(as.data.frame(ek), file.path(outdir, sprintf("%s_ORA_KEGG_%s.csv", prefix, dir)), row.names = FALSE)
          S(sprintf("   KEGG significativos: %d", nrow(as.data.frame(ek)))) }
      }
    }
  }
  if (gsea) {
    gl <- tt$t; names(gl) <- rownames(tt); gl <- sort(gl[is.finite(gl)], decreasing = TRUE)
    gg <- tryCatch(gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", minGSSize = 10, maxGSSize = 500,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE),
                   error = function(e) { S("   gseGO fallo: ", conditionMessage(e)); NULL })
    if (!is.null(gg)) {
      d <- as.data.frame(gg); write.csv(d, file.path(outdir, sprintf("%s_GSEA_GOBP.csv", prefix)), row.names = FALSE)
      S(sprintf("  [%s] GSEA GO BP (ranking = t moderado): %d vias (NES>0: %d; NES<0: %d)", prefix, nrow(d), sum(d$NES > 0), sum(d$NES < 0)))
      if (nrow(d)) { d <- d[order(-abs(d$NES)), ]; for (i in seq_len(min(8, nrow(d)))) S(sprintf("     NES=%+.2f q=%.2g  %s", d$NES[i], d$p.adjust[i], substr(d$Description[i], 1, 70))) }
      plot_top(gg, file.path(outdir, sprintf("%s_GSEA_GOBP.png", prefix)), sprintf("%s: GSEA GO BP", prefix))
    }
    if (kegg) {
      ent <- suppressWarnings(tryCatch(bitr(names(gl), "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db), error = function(e) NULL))
      if (!is.null(ent)) {
        ent <- ent[!duplicated(ent$SYMBOL) & !duplicated(ent$ENTREZID), ]
        gk <- gl[ent$SYMBOL]; names(gk) <- ent$ENTREZID; gk <- sort(gk, decreasing = TRUE)
        gkk <- tryCatch(gseKEGG(gk, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE),
                        error = function(e) { S("   gseKEGG no disponible: ", substr(conditionMessage(e), 1, 80)); NULL })
        if (!is.null(gkk)) { write.csv(as.data.frame(gkk), file.path(outdir, sprintf("%s_GSEA_KEGG.csv", prefix)), row.names = FALSE)
          S(sprintf("  [%s] GSEA KEGG: %d vias", prefix, nrow(as.data.frame(gkk)))) }
      }
    }
  }
  invisible(res)
}

# listas de marcadores (curadas a mano; revisar contra van Groningen 2017 / Boeva 2017 antes de citar)
MARKERS_MES  <- c("AXL","VIM","FN1","TGFBI","CAV1","CD44","SERPINE1","THY1","LGALS1","ANXA1","CAV2","S100A6","MGP","CCN1","CCN2",
                  "TAGLN","COL1A2","COL4A1","NT5E","PRRX1","YAP1","WWTR1","SNAI2","ZEB1")
MARKERS_ADRN <- c("PHOX2B","PHOX2A","HAND2","GATA3","ISL1","TH","DBH","DLK1","CHGA","HAND1","SYT1","STMN2","NEFL","GAP43",
                  "ELAVL4","ALK","CHRNA3","SLC18A2","PRPH","TFAP2B")
marker_summary <- function(tt, genes, label) {
  g <- genes[genes %in% rownames(tt)]
  S(sprintf("  %s: %d/%d genes presentes; FDR<0.05 con logFC>0: %d; con logFC<0: %d",
            label, length(g), length(genes), sum(tt[g, "adj.P.Val"] < FDR_CUT & tt[g, "logFC"] > 0), sum(tt[g, "adj.P.Val"] < FDR_CUT & tt[g, "logFC"] < 0)))
}
