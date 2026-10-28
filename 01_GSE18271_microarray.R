# =============================================================================
# 01_GSE18271_microarray.R   (Affymetrix HG-U133 Plus 2.0; 22 tumores: 12 ganglioneuroma (GN) y 10 ganglioneuroblastoma (GNB))
#   QC -> RMA -> anotacion y colapso sonda->gen -> limma (GNB vs GN) -> enriquecimiento
#   Incluye la reproduccion del analisis ORIGINAL (a nivel de sonda, sin filtro) para comparar con la memoria.
# Entradas : Datos/18271/GSE18271_series_matrix.txt  y  GSM4563xx.CEL.gz
# Salidas  : resultados/GSE18271/
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("limma", "AnnotationDbi", "ggplot2"))
if (!file.exists(file.path(RES, "GSE18271", "RMA_expr.rds"))) need(c("affy", "hgu133plus2.db"))
suppressPackageStartupMessages({ library(limma); library(AnnotationDbi) })
start_log("GSE18271")
D <- file.path(DATA, "18271")

# -----------------------------------------------------------------------------
SH("1. Metadatos (se leen de GEO, no se escriben a mano)")
meta <- parse_series_matrix(file.path(D, "GSE18271_series_matrix.txt"))
cat("Columnas:", paste(names(meta), collapse = ", "), "\n")
meta$group <- factor(ifelse(grepl("blastoma", meta$tissue), "GNB", "GN"), levels = c("GN", "GNB"))
rownames(meta) <- meta$geo_accession
S("Muestras: ", nrow(meta), " | ", paste(names(table(meta$group)), table(meta$group), collapse = " / "))
S("Orden de grupos (debe coincidir con el vector escrito a mano en el script antiguo): ", paste(tolower(as.character(meta$group)), collapse = ","))
old_group <- c("gn","gn","gnb","gn","gn","gn","gnb","gnb","gn","gnb","gnb","gn","gnb","gnb","gn","gn","gn","gnb","gnb","gn","gnb","gn")
S("Coincide con el vector antiguo: ", identical(tolower(as.character(meta$group)), old_group))

# -----------------------------------------------------------------------------
SH("2. Lectura de CEL y normalizacion RMA")
rds <- OUT("RMA_expr.rds")
if (file.exists(rds)) { expr <- readRDS(rds); S("(usando RMA cacheado: ", rds, ")") } else {
  suppressPackageStartupMessages(library(affy))
  files <- file.path(D, paste0(meta$geo_accession, ".CEL.gz")); stopifnot(all(file.exists(files)))
  raw <- ReadAffy(filenames = files, sampleNames = meta$geo_accession)
  S("Plataforma (cdf): ", cdfName(raw), " | arrays: ", ncol(exprs(raw)))
  png(OUT("QC_boxplot_crudo.png"), 1800, 1300, res = 200); boxplot(raw, las = 2, col = "lightblue", main = "Intensidades crudas (log2)"); dev.off()
  png(OUT("QC_densidad_crudo.png"), 1600, 1200, res = 200); hist(raw, main = "Densidades, datos crudos"); dev.off()
  eset <- rma(raw); expr <- exprs(eset); saveRDS(expr, rds)
}
write.csv(round(expr, 4), OUT("RMA_normalizado.csv"))
S("Matriz RMA: ", nrow(expr), " sondas x ", ncol(expr), " arrays")
png(OUT("QC_boxplot_RMA.png"), 1800, 1300, res = 200); boxplot(expr, las = 2, col = "lightgreen", main = "Tras RMA"); dev.off()
rle <- expr - apply(expr, 1, median)
png(OUT("QC_RLE.png"), 1800, 1300, res = 200); boxplot(rle, las = 2, col = "grey80", outline = FALSE, main = "RLE (log-expresion relativa a la mediana)"); abline(h = 0, col = "red"); dev.off()
cm <- cor(expr); mc <- (rowSums(cm) - 1) / (ncol(cm) - 1)
S("Correlacion media array-array: min ", round(min(mc), 4), " (", names(which.min(mc)), ") | max ", round(max(mc), 4),
  " | arrays a mas de 3 DE por debajo de la media: ", sum(mc < mean(mc) - 3 * sd(mc)))
S("Mediana del RLE por array (cercano a 0 es bueno): rango ", paste(round(range(apply(rle, 2, median)), 3), collapse = " a "))
ve <- pca_plot(expr, meta$group, OUT("PCA.png"), "PCA GSE18271 (todas las sondas)")
S("PCA: PC1 ", round(ve[1], 1), "%, PC2 ", round(ve[2], 1), "%")
pcs <- prcomp(t(expr), scale. = TRUE)$x; S("Separacion GN/GNB: PC1 p=", signif(t.test(pcs[, 1] ~ meta$group)$p.value, 3), " | PC2 p=", signif(t.test(pcs[, 2] ~ meta$group)$p.value, 3))

# -----------------------------------------------------------------------------
SH("3. Anotacion sonda -> simbolo")
probe2symbol <- function(p) {
  a <- suppressMessages(AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db, keys = p, keytype = "PROBEID", columns = c("PROBEID", "SYMBOL")))
  a <- a[!is.na(a$SYMBOL), ]; n <- table(a$PROBEID)
  a <- a[a$PROBEID %in% names(n)[n == 1], ]            # fuera sondas que mapean a >1 simbolo (ambiguas)
  setNames(a$SYMBOL, a$PROBEID)[p]
}
sym <- probe2symbol(rownames(expr)); names(sym) <- rownames(expr)
S("Sondas totales: ", nrow(expr), " | con un unico simbolo: ", sum(!is.na(sym)), " | sin simbolo o ambiguas: ", sum(is.na(sym)),
  " | simbolos unicos: ", length(unique(sym[!is.na(sym)])))

# -----------------------------------------------------------------------------
SH("4. ANALISIS ORIGINAL a nivel de SONDA (sin filtro) -> para comparar con la memoria")
des <- model.matrix(~ group, data = meta); colnames(des) <- c("Intercept", "GNB_vs_GN")
tt_p <- topTable(eBayes(lmFit(expr, des)), coef = "GNB_vs_GN", number = Inf); tt_p$Symbol <- sym[rownames(tt_p)]
dp <- tt_p[tt_p$adj.P.Val < FDR_CUT & abs(tt_p$logFC) > LFC_CUT, ]
S("Sondas con FDR<0.05 y |logFC|>1: ", nrow(dp), " (up ", sum(dp$logFC > 0), ", down ", sum(dp$logFC < 0), ") ; simbolos distintos entre ellas: ", length(unique(na.omit(dp$Symbol))))
S("   memoria: 2076 'genes' DEG (comprobar si eran sondas o genes)")
oldf <- file.path(D, "DEG_results_GSE18271.csv")
if (file.exists(oldf)) { o <- read.csv(oldf, row.names = 1); cp <- intersect(rownames(o), rownames(tt_p))
  S("Comparacion con tu DEG_results_GSE18271.csv: sondas comunes ", length(cp), "; correlacion logFC = ", round(cor(o[cp, "logFC"], tt_p[cp, "logFC"]), 6),
    "; DEGs en tu CSV = ", sum(o$adj.P.Val < FDR_CUT & abs(o$logFC) > LFC_CUT)) }
S("Sondas de los genes candidatos (nivel sonda):"); pc <- tt_p[tt_p$Symbol %in% CANDIDATES & !is.na(tt_p$Symbol), c("Symbol", "logFC", "AveExpr", "adj.P.Val")]; STAB(pc[order(pc$Symbol, pc$adj.P.Val), ], 3)

# -----------------------------------------------------------------------------
SH("5. ANALISIS CORREGIDO a nivel de GEN")
keepA <- !grepl("^AFFX", rownames(expr)) & !is.na(sym)
cutoff <- quantile(expr, 0.25); nmin <- min(table(meta$group))
keepE <- rowSums(expr > cutoff) >= nmin
S("Filtro: sonda con simbolo, no control, y expresion > percentil 25 global (", round(cutoff, 2), " log2) en >= ", nmin, " arrays")
ex <- expr[keepA & keepE, ]; sy <- sym[rownames(ex)]
S("Sondas que pasan: ", nrow(ex))
# colapso: por cada simbolo se queda la sonda de mayor expresion media
best <- tapply(seq_len(nrow(ex)), sy, function(i) i[which.max(rowMeans(ex[i, , drop = FALSE]))])
Ge <- ex[unlist(best), ]; rownames(Ge) <- names(best)
S("Genes tras colapsar (sonda de mayor expresion media por gen): ", nrow(Ge))
fit <- eBayes(lmFit(Ge, des), robust = TRUE)
tt <- topTable(fit, coef = "GNB_vs_GN", number = Inf); write.csv(tt, OUT("DE_GNB_vs_GN_gen.csv"))
S("Genes testados ", nrow(tt)); STAB(t(as.data.frame(deg_counts(tt))), 7)
S("Candidatos (logFC = GNB vs GN; positivo = mas alto en GNB):"); STAB(cand_table(tt), 3)
S("Memoria: ARPP21 +2.97, DLG2 +2.72")
volcano_plot(tt, "GNB vs GN (GSE18271, nivel gen)", OUT("volcano.png"), xlab = "log2FC (GNB / GN)")
heatmap_png(Ge[rownames(tt)[1:50], ], data.frame(Grupo = meta$group, row.names = colnames(Ge)), OUT("heatmap_top50.png"), "Top 50 genes (GNB vs GN)")
heatmap_png(cm, data.frame(Grupo = meta$group, row.names = colnames(cm)), OUT("correlacion_arrays.png"), "Correlacion entre arrays", scale = "none")
write.csv(tt[1:20, c("logFC", "AveExpr", "adj.P.Val")], OUT("top20.csv"))
S("Top 10 genes por FDR:"); STAB(tt[1:10, c("logFC", "AveExpr", "adj.P.Val")], 3)

SH("6. Posible confusion con la proporcion de estroma schwanniano (marcadores)")
SCHWANN <- c("S100B", "SOX10", "PLP1", "MPZ", "GFAP", "CDH19", "ERBB3", "PMP22", "MBP", "NGFR", "SOX2")
NEURONAL <- c("PHOX2B", "STMN2", "NEFL", "SYP", "GAP43", "ELAVL4", "TUBB3", "DBH", "TH", "HAND2", "GATA3")
sh <- function(g, lab) { g <- g[g %in% rownames(tt)]; d <- tt[g, c("logFC", "adj.P.Val")]; S(lab, ": mediana logFC = ", round(median(d$logFC), 2), "; ", sum(d$adj.P.Val < FDR_CUT & d$logFC > 0), " suben y ", sum(d$adj.P.Val < FDR_CUT & d$logFC < 0), " bajan (FDR<0.05) de ", length(g)); d }
d1 <- sh(SCHWANN, "Marcadores de celula de Schwann"); STAB(d1, 3); d2 <- sh(NEURONAL, "Marcadores neuronales/adrenergicos"); STAB(d2, 3)
S("Lectura: logFC negativo en marcadores de Schwann = GN tiene MAS estroma schwanniano que GNB; entonces la firma GNB vs GN refleja en parte composicion celular (neuroblastos/ganglionares vs estroma), no solo 'diferenciacion' intrinseca.")

SH("7. Enriquecimiento funcional (GNB vs GN)")
run_enrichment(tt, "GNBvsGN")
end_log()
