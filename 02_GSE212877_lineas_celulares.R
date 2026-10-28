# =============================================================================
# 02_GSE212877_lineas_celulares.R   (RNA-seq de 8 lineas celulares de neuroblastoma, controles sin tratar)
# Contrastes que se ejecutan (todos con voom + limma):
#   A  TFM ORIGINAL : SHEP (MES) vs 4 lineas ADRN (55N, BE2, SY5Y, WN)        -> reproduce la memoria
#   A0 (version antigua sin voom, solo para comparar)
#   B  NULOS        : cada linea ADRN contra las otras tres ADRN (mismo fenotipo)  -> cuanto "DEG" produce la identidad de linea
#   C  PAREADO      : tipo S vs tipo N dentro de pares de lineas (55N/5S, WN/WS, SY5Y/SHEP), bloque = par
#   D  PAREADO+mixto: igual que C pero con linea como efecto aleatorio (duplicateCorrelation) -> mas conservador
#   E  S (SHEP,5S,WS) vs N (resto, 5 lineas) sin parear
#   F  Robustez: C quitando cada par por turnos
# Entradas : Datos/212877/GSE212877_series_matrix.txt y los GSM*_C[1-3]_RNASeq.txt.gz
# Salidas  : resultados/GSE212877/
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("limma", "ggplot2"))
suppressPackageStartupMessages(library(limma))
start_log("GSE212877")
D <- file.path(DATA, "212877")

# -----------------------------------------------------------------------------
SH("1. Metadatos y lectura de conteos")
meta <- parse_series_matrix(file.path(D, "GSE212877_series_matrix.txt"))
cat("Columnas:", paste(names(meta), collapse = ", "), "\n")
S("Serie con ", nrow(meta), " muestras; tratamientos: ", paste(names(table(meta$treatment)), table(meta$treatment), collapse = " | "))
ctrl <- meta[meta$treatment == "none", ]
stopifnot(nrow(ctrl) == 24)
ctrl$line <- sub(" .*", "", ctrl$title)
S("Controles sin tratar: ", nrow(ctrl), " (", paste(names(table(ctrl$line)), table(ctrl$line), collapse = ", "), ")")

# Informacion de cada linea (fuente: nombres/clasificacion de GEO + literatura; REVISAR antes de citar)
line_info <- data.frame(
  line  = c("55N", "5S", "BE2", "NBLS", "SHEP", "SY5Y", "WN", "WS"),
  geo   = c("LA1-55n", "LA1-5s", "SK-N-BE2", "NBL-S", "SHEP", "SH-SY5Y", "NBL-W-N", "NBL-W-S"),
  type  = c("N", "S", "N", "N", "S", "N", "N", "S"),            # N = neuroblastico/adrenergico ; S = sustrato-adherente/mesenquimal-like
  pair  = c("LA1", "LA1", NA, NA, "SKNSH", "SKNSH", "NBLW", "NBLW"),  # LA1-55n/5s ; SK-N-SH -> SY5Y/SHEP ; NBL-W-N/S
  mycn  = c("amp", "amp", "amp", "overexp", "nonamp", "nonamp", "amp", "amp"),
  tfm   = c("ADRN", "hibrido", "ADRN", "excluida", "MES", "ADRN", "ADRN", "hibrido"))  # asignacion usada en la memoria
STAB(line_info)
chk <- tapply(ctrl$source, ctrl$line, function(x) unique(x)); S("Comprobacion nombre GEO por linea: ", paste(names(chk), unlist(chk), sep = "=", collapse = "; "))
stopifnot(all(unlist(chk)[line_info$line] == line_info$geo))

read_counts <- function(geo) {
  f <- list.files(D, pattern = paste0("^", geo, "_.*RNASeq\\.txt\\.gz$"), full.names = TRUE)
  stopifnot(length(f) == 1)
  x <- read.table(gzfile(f), header = FALSE, col.names = c("gene", "count"), stringsAsFactors = FALSE, quote = "", comment.char = "")
  setNames(x$count, x$gene)
}
cl <- lapply(ctrl$geo_accession, read_counts)
stopifnot(length(unique(lapply(cl, names))) == 1)
counts <- do.call(cbind, cl); colnames(counts) <- paste0(ctrl$line, "_", sub(".*rep", "r", ctrl$title))
rownames(counts) <- names(cl[[1]])
S("Matriz de conteos: ", nrow(counts), " genes x ", ncol(counts), " muestras; genes duplicados: ", sum(duplicated(rownames(counts))))
S("Lecturas por muestra (millones): ", paste(round(range(colSums(counts)) / 1e6, 1), collapse = " - "))
line <- factor(ctrl$line, levels = line_info$line); names(line) <- colnames(counts)
type <- factor(line_info$type[match(line, line_info$line)], levels = c("N", "S"))
pair <- line_info$pair[match(line, line_info$line)]

# -----------------------------------------------------------------------------
SH("2. Filtrado, normalizacion TMM y exploracion")
tmm_factors <- function(cn) {
  if (has("edgeR")) return(edgeR::calcNormFactors(edgeR::DGEList(cn))$samples$norm.factors)
  lib <- colSums(cn); ref <- which.min(abs(lib - mean(lib)))
  nf <- sapply(seq_len(ncol(cn)), function(i) { if (i == ref) return(1)
    o <- cn[, i]; r <- cn[, ref]; no <- sum(o); nr <- sum(r); ok <- o > 0 & r > 0
    M <- log2((o[ok] / no) / (r[ok] / nr)); A <- 0.5 * log2((o[ok] / no) * (r[ok] / nr))
    w <- 1 / ((no - o[ok]) / no / o[ok] + (nr - r[ok]) / nr / r[ok])
    qM <- quantile(M, c(.3, .7)); qA <- quantile(A, c(.05, .95)); s <- M >= qM[1] & M <= qM[2] & A >= qA[1] & A <= qA[2]
    2^(sum(w[s] * M[s]) / sum(w[s])) })
  nf / exp(mean(log(nf)))
}
if (has("edgeR")) {
  S("Usando edgeR (filterByExpr + TMM), igual que el analisis original")
  dge <- edgeR::DGEList(counts, group = line); keep <- edgeR::filterByExpr(dge, group = line)
} else {
  S("AVISO: edgeR no instalado; uso un filtro y un TMM propios (los numeros pueden diferir ligeramente de la memoria)")
  cpm_ <- t(t(counts) / colSums(counts) * 1e6); keep <- rowSums(cpm_ > 10 / median(colSums(counts)) * 1e6) >= 3 & rowSums(counts) >= 15
}
cnt <- counts[keep, ]
S("Genes tras filtro de expresion: ", nrow(cnt), " de ", nrow(counts))
nf <- tmm_factors(cnt); libeff <- colSums(cnt) * nf
logcpm <- log2(t(t(cnt + 0.5) / (libeff + 1) * 1e6))
write.csv(data.frame(muestra = colnames(cnt), linea = line, lib = colSums(cnt), norm_factor = nf), OUT("muestras_QC.csv"), row.names = FALSE)

ve <- pca_plot(logcpm, line, OUT("PCA_8_lineas.png"), "PCA - 8 lineas celulares (controles)")
S("PCA: PC1 ", round(ve[1], 1), "% , PC2 ", round(ve[2], 1), "%")
png(OUT("MDS_8_lineas.png"), 1500, 1100, res = 200); plotMDS(logcpm, labels = line, col = as.integer(line)); dev.off()
cor_mat <- cor(logcpm)
heatmap_png(cor_mat, data.frame(linea = line, tipo = type, row.names = colnames(cnt)), OUT("correlacion_muestras.png"), "Correlacion entre muestras", scale = "none")
png(OUT("dendrograma.png"), 1500, 1000, res = 200); plot(hclust(dist(t(logcpm))), main = "Clustering jerarquico (logCPM)"); dev.off()
S("Correlacion media entre replicas de una misma linea: ", round(mean(sapply(levels(line), function(l) { m <- cor_mat[line == l, line == l]; mean(m[upper.tri(m)]) })), 4),
  " | entre lineas distintas: ", round(mean(cor_mat[outer(line, line, "!=")]), 4))

# -----------------------------------------------------------------------------
vcontrast <- function(sel, des, coef, block = NULL) {
  cn <- cnt[, sel, drop = FALSE]; lib <- colSums(cn) * tmm_factors(cn)
  v <- voom(cn, des, lib.size = lib)
  if (!is.null(block)) {
    dc <- duplicateCorrelation(v, des, block = block)
    v <- voom(cn, des, lib.size = lib, block = block, correlation = dc$consensus)
    dc <- duplicateCorrelation(v, des, block = block)
    fit <- lmFit(v, des, block = block, correlation = dc$consensus)
    attr(fit, "dupcor") <- dc$consensus
  } else fit <- lmFit(v, des)
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = coef, number = Inf)
  attr(tt, "dupcor") <- attr(fit, "dupcor"); tt
}
summ <- function(tt, label) { S(sprintf("%-46s", label), paste(names(deg_counts(tt)), deg_counts(tt), collapse = "  ")) }
RES_ALL <- list()

SH("3. Contraste A: TFM ORIGINAL, SHEP (MES) vs ADRN (55N, BE2, SY5Y, WN)")
ph <- factor(line_info$tfm[match(line, line_info$line)], levels = c("ADRN", "MES", "hibrido", "excluida"))
selA <- ph %in% c("ADRN", "MES"); phA <- droplevels(ph[selA]); S("Muestras: ", paste(names(table(phA)), table(phA), collapse = " / "))
desA <- model.matrix(~ phA)
ttA <- vcontrast(selA, desA, "phAMES"); RES_ALL$A_TFM_SHEP_vs_ADRN4 <- ttA
write.csv(ttA, OUT("A_SHEPvsADRN4_voom.csv")); summ(ttA, "A  SHEP vs ADRN4 (voom)")
S("   memoria: 6005 DEGs (FDR<0.05 & |logFC|>1)"); STAB(cand_table(ttA), 3)
marker_summary(ttA, MARKERS_MES, "marcadores MES"); marker_summary(ttA, MARKERS_ADRN, "marcadores ADRN")
# A0: sin voom (analisis antiguo)
fit0 <- eBayes(lmFit(logcpm[, selA], desA)); tt0 <- topTable(fit0, coef = "phAMES", number = Inf); summ(tt0, "A0 SHEP vs ADRN4 (logCPM sin voom, antiguo)")
# comparacion con los CSV antiguos del usuario (si estan)
old <- file.path(D, "DEGs_ADRN_vs_MES_voom.csv")
if (file.exists(old)) { o <- read.csv(old, row.names = 1); cg <- intersect(rownames(o), rownames(ttA))
  S("Comparacion con tu DEGs_ADRN_vs_MES_voom.csv: genes comunes ", length(cg), "; correlacion de logFC = ", round(cor(o[cg, "logFC"], ttA[cg, "logFC"]), 4),
    "; DEGs en tu CSV = ", sum(o$adj.P.Val < FDR_CUT & abs(o$logFC) > LFC_CUT)) }
volcano_plot(ttA, "SHEP (MES) vs ADRN4 (voom)", OUT("volcano_A_SHEP_vs_ADRN4.png"), xlab = "log2FC (MES / ADRN)")

SH("4. Contraste B: NULOS (una linea ADRN vs las otras ADRN)")
ADRN4 <- c("55N", "BE2", "SY5Y", "WN"); nulls <- NULL
for (L in ADRN4) { sel <- line %in% ADRN4; g <- factor(ifelse(line[sel] == L, "A", "B"), levels = c("B", "A"))
  tt <- vcontrast(sel, model.matrix(~ g), "gA"); dc <- deg_counts(tt)
  nulls <- rbind(nulls, data.frame(linea = L, FDR_y_LFC = unname(dc["FDR_y_LFC"]), up = unname(dc["up"]), down = unname(dc["down"]))) }
STAB(nulls, 7); S("Referencia: contraste A (SHEP vs ADRN4) = ", deg_counts(ttA)["FDR_y_LFC"], " DEGs. Si los nulos dan una fraccion grande de esa cifra, buena parte de A es identidad de linea.")
write.csv(nulls, OUT("B_nulos.csv"), row.names = FALSE)

SH("5. Contraste C: PAREADO S vs N (bloque = par), 3 pares")
pairs3 <- c("55N", "5S", "WN", "WS", "SY5Y", "SHEP"); selC <- line %in% pairs3
pairC <- factor(pair[selC]); typeC <- droplevels(type[selC])
S("Pares: LA1 (55N/5S), NBLW (WN/WS), SKNSH (SY5Y/SHEP); replicas por linea: 3")
desC <- model.matrix(~ pairC + typeC)
ttC <- vcontrast(selC, desC, "typeCS"); RES_ALL$C_paired <- ttC; write.csv(ttC, OUT("C_pareado_S_vs_N.csv")); summ(ttC, "C  pareado S vs N")
STAB(cand_table(ttC), 3); marker_summary(ttC, MARKERS_MES, "marcadores MES"); marker_summary(ttC, MARKERS_ADRN, "marcadores ADRN")
volcano_plot(ttC, "S vs N pareado (voom)", OUT("volcano_C_pareado.png"), xlab = "log2FC (S / N)")

SH("6. Contraste D: pareado + linea como efecto aleatorio (duplicateCorrelation)")
blockD <- factor(as.character(line[selC]))
ttD <- vcontrast(selC, desC, "typeCS", block = blockD); RES_ALL$D_paired_mixed <- ttD; write.csv(ttD, OUT("D_pareado_mixto.csv"))
S("Correlacion intra-linea estimada (consensus): ", round(attr(ttD, "dupcor"), 3)); summ(ttD, "D  pareado + efecto aleatorio de linea")
STAB(cand_table(ttD), 3); marker_summary(ttD, MARKERS_MES, "marcadores MES"); marker_summary(ttD, MARKERS_ADRN, "marcadores ADRN")
volcano_plot(ttD, "S vs N pareado, efecto aleatorio de linea", OUT("volcano_D_pareado_mixto.png"), xlab = "log2FC (S / N)")
cg <- intersect(rownames(ttA), rownames(ttD))
S("Concordancia A vs D: correlacion logFC = ", round(cor(ttA[cg, "logFC"], ttD[cg, "logFC"]), 3),
  " | DEGs de A tambien DEG en D con el mismo signo: ", sum(ttA[cg, "adj.P.Val"] < FDR_CUT & abs(ttA[cg, "logFC"]) > LFC_CUT & ttD[cg, "adj.P.Val"] < FDR_CUT & abs(ttD[cg, "logFC"]) > LFC_CUT & sign(ttA[cg, "logFC"]) == sign(ttD[cg, "logFC"])))

SH("7. Contraste E: S (SHEP, 5S, WS) vs N (55N, BE2, NBLS, SY5Y, WN), sin parear")
typeE <- factor(ifelse(line %in% c("SHEP", "5S", "WS"), "S", "N"), levels = c("N", "S"))
ttE <- vcontrast(rep(TRUE, ncol(cnt)), model.matrix(~ typeE), "typeES"); write.csv(ttE, OUT("E_S_vs_N_sin_parear.csv")); summ(ttE, "E  S vs N sin parear"); STAB(cand_table(ttE), 3)

SH("8. Contraste F: robustez, pareado quitando cada par")
Fres <- NULL
for (drop in c("LA1", "NBLW", "SKNSH")) {
  sel <- line %in% pairs3 & !(pair %in% drop); pp <- factor(pair[sel]); tp <- droplevels(type[sel])
  tt <- vcontrast(sel, model.matrix(~ pp + tp), "tpS"); ct <- cand_table(tt)
  Fres <- rbind(Fres, data.frame(quitado = drop, gen = CANDIDATES, logFC = ct$logFC, FDR = ct$adj.P.Val))
}
STAB(Fres, 3); write.csv(Fres, OUT("F_robustez_quitando_pares.csv"), row.names = FALSE)

SH("9. Detectabilidad de los candidatos (cuanto se expresan realmente en las lineas)")
S("Conteos crudos medios por linea (", round(mean(colSums(counts)) / 1e6, 1), " M lecturas/muestra):")
STAB(round(sapply(levels(line), function(l) rowMeans(counts[CANDIDATES, line == l, drop = FALSE])), 1), 1)
S("CPM medios por linea:")
cpm_all <- t(t(counts) / colSums(counts) * 1e6)
STAB(round(sapply(levels(line), function(l) rowMeans(cpm_all[CANDIDATES, line == l, drop = FALSE])), 2), 2)
S("Genes con menos de 10 lecturas en todas las muestras no son interpretables (ARPP21 es el caso a mirar).")

SH("10. Resumen de contrastes")
tab <- t(sapply(RES_ALL, deg_counts)); STAB(tab, 7)
S("Candidatos en cada contraste (logFC | FDR):")
for (nm in names(RES_ALL)) { ct <- cand_table(RES_ALL[[nm]]); S(sprintf("  %-22s %s", nm, paste(sprintf("%s %+.2f|%.2g", CANDIDATES, ct$logFC, ct$adj.P.Val), collapse = "   "))) }

SH("11. Enriquecimiento funcional")
S("-- Contraste A (el de la memoria) --"); run_enrichment(ttA, "A_SHEPvsADRN4")
S("-- Contraste D (pareado, conservador) --"); run_enrichment(ttD, "D_pareado_mixto")
SH("12. Sensibilidad del contraste D (referencia): estabilidad del modelo y genes de conteos muy bajos")
cmp_de <- function(a, b, nm) { g <- intersect(rownames(a), rownames(b)); sa <- a[g, "adj.P.Val"] < FDR_CUT & abs(a[g, "logFC"]) > LFC_CUT; sb <- b[g, "adj.P.Val"] < FDR_CUT & abs(b[g, "logFC"]) > LFC_CUT
  S(sprintf("  D vs %s: correlacion logFC = %.3f | Spearman de t = %.3f | DEG de D tambien DEG en %s con el mismo signo: %d de %d",
            nm, cor(a[g, "logFC"], b[g, "logFC"]), cor(a[g, "t"], b[g, "t"], method = "spearman"), nm, sum(sa & sb & sign(a[g, "logFC"]) == sign(b[g, "logFC"])), sum(sa))) }
S("(a) Alternativa con la linea como efecto fijo (C: bloque por par) y sin parear (E)")
cmp_de(ttD, ttC, "C"); cmp_de(ttD, ttE, "E")
S("(b) Genes con |log2FC| > 8 en D (conteos muy bajos en uno de los grupos)")
ext <- abs(ttD$logFC) > 8
S(sprintf("  %d genes con |log2FC| > 8 (%d DEG; %d con expresion media < 1 log2 CPM; mas en S: %d, mas en N: %d)", sum(ext), sum(ext & ttD$adj.P.Val < FDR_CUT), sum(ext & ttD$AveExpr < 1), sum(ext & ttD$logFC > 0), sum(ext & ttD$logFC < 0)))
if (all(vapply(c("clusterProfiler", "org.Hs.eg.db"), has, logical(1)))) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  uni <- rownames(ttD)
  for (dir in c("up", "down")) {
    sel <- ttD$adj.P.Val < FDR_CUT & (if (dir == "up") ttD$logFC > LFC_CUT else ttD$logFC < -LFC_CUT)
    g_all <- rownames(ttD)[sel]; g_ex <- rownames(ttD)[sel & !ext]
    e1 <- as.data.frame(enrichGO(g_all, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = uni, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
    e2 <- as.data.frame(enrichGO(g_ex, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = uni, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
    top <- head(e1$ID, 20)
    S(sprintf("  ORA %s: %d genes (%d sin los extremos) | terminos GO BP q<0,05: %d con todos, %d sin los extremos, %d en comun | de los 20 primeros, siguen significativos sin los extremos: %d",
              if (dir == "up") "mas en S" else "mas en N", length(g_all), length(g_ex), nrow(e1), nrow(e2), length(intersect(e1$ID, e2$ID)), sum(top %in% e2$ID)))
    S("     primeros terminos sin los extremos: ", paste(head(e2$Description, 4), collapse = "; "))
  }
} else S("  (ORA omitido: falta clusterProfiler/org.Hs.eg.db)")
end_log()
