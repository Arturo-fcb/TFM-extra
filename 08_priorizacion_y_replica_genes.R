# =============================================================================
# 08_priorizacion_y_replica_genes.R   SELECCION DE CANDIDATOS TRAS LOS RESULTADOS GLOBALES + REPLICA EXTERNA
#   Pregunta: ?que genes se asocian con la supervivencia global DESPUES de ajustar por edad, MYCN y estadio
#   en SEQC (descubrimiento, GSE62564) y se replican, con el mismo sentido, en TARGET-NBL (validacion)?
#   - Criterios fijados ANTES de mirar ningun gen concreto (ver CRITERIOS). Ningun gen se elige a mano.
#   - Los cuatro genes de la version inicial del TFM (ARPP21, DLG2, AXL, MYCN) se evaluan al final como caso
#     de estudio a posteriori; no se usan para definir nada.
# Entradas : resultados/cache_G_GSE62564.rds (script 04), Datos/62564 y Datos/49711 (series matrix),
#            resultados/TARGET_NBL/cache_TARGET.rds (script 06, version con cache)
# Salidas  : resultados/replica_genes/  (tablas CSV, fig_replica_genes.png, RESUMEN_replica_genes.txt)
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("survival", "ggplot2", "ggrepel"))
suppressPackageStartupMessages({ library(survival); library(ggplot2); library(ggrepel) })
start_log("replica_genes")
num <- function(x) suppressWarnings(as.numeric(x))
PT <- function(x) for (l in capture.output(print(x))) S(l)
D <- file.path(DATA, "62564")

CRITERIOS <- c("Modelo: Cox de supervivencia global ~ expresion del gen (z-score) + log(edad) + MYCN amplificado + estadio 4",
               "Descubrimiento: SEQC, FDR (BH, todos los genes) < 0,05",
               "Replica: TARGET-NBL, mismo sentido del efecto y p < 0,05 (ademas: FDR de BH dentro de los genes descubiertos)")
SH("0. Criterios de seleccion (fijados a priori)"); for (c in CRITERIOS) S("  - ", c)

# ---- 1. SEQC -----------------------------------------------------------------
SH("1. Datos SEQC (descubrimiento) y TARGET-NBL (validacion)")
cache <- file.path(RES, "cache_G_GSE62564.rds"); if (!file.exists(cache)) stop("Falta ", cache, " (script 04)")
G <- readRDS(cache)
meta <- parse_series_matrix(file.path(D, "GSE62564_series_matrix.txt")); meta$id <- sub(" \\[.*\\]$", "", meta$title)
cand <- c(file.path(DATA, "49711", "GSE49711_series_matrix.txt"), file.path(DATA, "49711", "GSE49711_series_matrix.txt.gz"),
          file.path(DATA, "GSE49711_series_matrix.txt"), file.path(DATA, "GSE49711_series_matrix.txt.gz"))
f49 <- cand[file.exists(cand)][1]; if (is.na(f49)) stop("No encuentro GSE49711_series_matrix.txt(.gz)")
m49 <- parse_series_matrix(f49); m49$id <- sub(" \\[.*\\]$", "", m49$title); m49 <- m49[match(meta$id, m49$id), ]
sa <- data.frame(id = meta$id, os_days = num(meta$os.day), os_event = num(meta$os.bin), efs_days = num(meta$efs.day), efs_event = num(meta$efs.bin), age_days = num(meta$age))
sa$mycn_amp <- ifelse(m49$mycn.status %in% c("0", "1"), as.integer(m49$mycn.status), NA)
sa$stage4 <- as.integer(trimws(m49$inss.stage) == "4"); sa$stage4[is.na(m49$inss.stage) | !nzchar(trimws(m49$inss.stage))] <- NA
sa <- sa[complete.cases(sa) & sa$os_days > 0, ]; rownames(sa) <- sa$id; sa$logage <- log1p(sa$age_days)
Gs <- G[rowSums(G >= 0) >= 0.2 * ncol(G), sa$id]
tf <- file.path(RES, "TARGET_NBL", "cache_TARGET.rds"); if (!file.exists(tf)) stop("Falta ", tf, ": vuelve a ejecutar el script 06 (version con cache)")
T0 <- readRDS(tf); st <- T0$sv
st <- st[complete.cases(st[, c("os_days", "os_event", "efs_days", "efs_event", "age_days", "mycn_amp", "stage4")]) & st$os_days > 0, ]
st$logage <- log1p(st$age_days); Gt <- T0$Gf[, st$case]
S("SEQC: ", nrow(sa), " pacientes, ", sum(sa$os_event), " muertes, ", nrow(Gs), " genes | TARGET-NBL: ", nrow(st), " pacientes, ", sum(st$os_event), " muertes, ", nrow(Gt), " genes")
genes <- intersect(rownames(Gs), rownames(Gt)); S("Genes comunes: ", length(genes))

# ---- 2. Cribado Cox gen a gen ------------------------------------------------
scan <- function(Gm, d, time, ev) {
  Z <- t(scale(t(Gm[genes, d[[1]], drop = FALSE])))
  S0 <- Surv(d[[time]], d[[ev]]); out <- matrix(NA_real_, length(genes), 3, dimnames = list(genes, c("coef", "z", "p")))
  x1 <- d$logage; x2 <- d$mycn_amp; x3 <- d$stage4
  for (i in seq_along(genes)) {
    m <- tryCatch(coxph(S0 ~ Z[i, ] + x1 + x2 + x3), error = function(e) NULL, warning = function(w) suppressWarnings(coxph(S0 ~ Z[i, ] + x1 + x2 + x3)))
    if (!is.null(m)) { s <- summary(m)$coefficients[1, ]; out[i, ] <- c(s["coef"], s["z"], s["Pr(>|z|)"]) }
  }
  out <- as.data.frame(out); out$fdr <- p.adjust(out$p, "BH"); out
}
SH("2. Cribado Cox de todos los genes comunes (HR por 1 DE de expresion; ajustado por edad, MYCN y estadio 4)")
t0 <- Sys.time()
os_s <- scan(Gs, sa, "os_days", "os_event"); S("  SEQC OS listo")
os_t <- scan(Gt, st, "os_days", "os_event"); S("  TARGET OS listo")
ef_s <- scan(Gs, sa, "efs_days", "efs_event"); S("  SEQC EFS listo")
ef_t <- scan(Gt, st, "efs_days", "efs_event"); S("  TARGET EFS listo | tiempo: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# ---- 3. Replica --------------------------------------------------------------
repl <- function(a, b, etiqueta) {
  ok <- is.finite(a$z) & is.finite(b$z); rho <- cor(a$z[ok], b$z[ok], method = "spearman")
  sig <- ok & a$fdr < FDR_CUT; same <- sig & sign(a$coef) == sign(b$coef); rep <- same & b$p < 0.05
  nsig <- sum(sig); frac <- sum(rep) / max(nsig, 1)
  bt <- if (nsig > 0) binom.test(sum(rep), nsig, p = 0.025, alternative = "greater")$p.value else NA   # esperado por azar: 2,5 % (mismo sentido y p < 0,05)
  S(sprintf("%s: genes con FDR<0,05 en SEQC = %d | mismo sentido en TARGET = %d (%.0f%%) | ademas p<0,05 en TARGET = %d (%.1f%%; esperado por azar 2,5%%; binomial p = %.2g) | Spearman de z entre cohortes (todos los genes) = %.3f",
            etiqueta, nsig, sum(same), 100 * sum(same) / max(nsig, 1), sum(rep), 100 * frac, bt, rho))
  data.frame(gen = rownames(a), coef_SEQC = a$coef, z_SEQC = a$z, p_SEQC = a$p, fdr_SEQC = a$fdr, coef_TARGET = b$coef, z_TARGET = b$z, p_TARGET = b$p,
             descubierto = sig, mismo_sentido = same, replicado = rep)
}
SH("3. Replica en TARGET-NBL")
R_os <- repl(os_s, os_t, "OS "); R_ef <- repl(ef_s, ef_t, "EFS")
R_os$fdr_TARGET_en_descubiertos <- NA_real_; i <- R_os$descubierto; R_os$fdr_TARGET_en_descubiertos[i] <- p.adjust(R_os$p_TARGET[i], "BH")
R_os$replicado_BH <- R_os$replicado & !is.na(R_os$fdr_TARGET_en_descubiertos) & R_os$fdr_TARGET_en_descubiertos < 0.05
S("Genes replicados en OS: p<0,05 en TARGET = ", sum(R_os$replicado), " | con FDR(BH) < 0,05 dentro de los descubiertos = ", sum(R_os$replicado_BH))
S("Genes replicados tambien en EFS (descubiertos y p<0,05 en TARGET, mismo sentido en EFS): ", sum(R_ef$replicado))
R_os$replicado_EFS <- R_ef$replicado[match(R_os$gen, R_ef$gen)]
write.csv(R_os, OUT("tabla_replica_genes_OS.csv"), row.names = FALSE); write.csv(R_ef, OUT("tabla_replica_genes_EFS.csv"), row.names = FALSE)
top <- R_os[R_os$replicado, ]; top <- top[order(top$p_TARGET), ]
S("Genes replicados (OS), ordenados por p en TARGET (primeros 40):")
PT(head(data.frame(gen = top$gen, z_SEQC = round(top$z_SEQC, 2), fdr_SEQC = signif(top$fdr_SEQC, 2), z_TARGET = round(top$z_TARGET, 2), p_TARGET = signif(top$p_TARGET, 2),
                   fdr_TARGET = signif(top$fdr_TARGET_en_descubiertos, 2), EFS = top$replicado_EFS), 40))
S("Sentido: de los replicados, ", sum(top$coef_SEQC > 0), " son de riesgo (HR>1) y ", sum(top$coef_SEQC < 0), " de proteccion (HR<1)")

# ---- 4. Enriquecimiento de los genes replicados ------------------------------
SH("4. Enriquecimiento funcional de los genes replicados (universo = genes comunes)")
if (all(vapply(c("clusterProfiler", "org.Hs.eg.db"), has, logical(1))) && nrow(top) >= 10) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  for (nm in c("riesgo", "proteccion")) {
    g <- top$gen[if (nm == "riesgo") top$coef_SEQC > 0 else top$coef_SEQC < 0]; S("  ", nm, ": ", length(g), " genes")
    if (length(g) < 10) next
    ego <- tryCatch(enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = genes, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05), error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego))) {
      ego2 <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) ego)
      write.csv(as.data.frame(ego2), OUT(sprintf("replicados_%s_ORA_GOBP_simplificado.csv", nm)), row.names = FALSE)
      d <- head(as.data.frame(ego2)[, c("Description", "p.adjust")], 8); for (j in seq_len(nrow(d))) S(sprintf("     %-70s q=%.2g", substr(d$Description[j], 1, 70), d$p.adjust[j]))
    } else S("   sin terminos significativos")
  }
} else S("  (omitido: falta clusterProfiler/org.Hs.eg.db o menos de 10 genes replicados)")

# ---- 4b. Solapamiento de los genes replicados con los programas y con la amplificacion de MYCN ---
SH("4b. Solapamiento de los genes de riesgo replicados con los programas biologicos y con la firma de amplificacion del oncogen MYCN")
fp <- file.path(RES, "GSE62564_score", "genes_programas.csv"); fm <- file.path(RES, "GSE62564_ajuste", "DE_MYCNamp_vs_noamp.csv")
if (file.exists(fp) && file.exists(fm)) {
  GP <- read.csv(fp, stringsAsFactors = FALSE); DM <- read.csv(fm, stringsAsFactors = FALSE, check.names = FALSE); DM$gen <- DM[[1]]
  rr <- R_os$gen[R_os$replicado & R_os$coef_SEQC > 0]
  conj <- c(split(GP$gen, GP$programa), list("MYCN_amp_up" = DM$gen[DM$adj.P.Val < 0.05 & DM$logFC > 1]))
  ov <- do.call(rbind, lapply(names(conj), function(k) {
    g <- intersect(conj[[k]], genes); a <- length(intersect(rr, g)); b <- length(rr) - a; c <- length(g) - a; d <- length(genes) - a - b - c
    data.frame(conjunto = k, genes_en_conjunto = length(g), solapan = a, pct_en_replicados = 100 * a / length(rr), pct_en_resto = 100 * (length(g) - a) / (length(genes) - length(rr)),
               OR = (a * d) / (b * c), p_Fisher = fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")$p.value)
  }))
  ov$fdr <- p.adjust(ov$p_Fisher, "BH"); write.csv(ov, OUT("solapamiento_replicados_riesgo.csv"), row.names = FALSE)
  S("  genes de riesgo replicados: ", length(rr), " (universo: ", length(genes), " genes comunes)")
  PT(data.frame(ov[, 1:3], pct_rep = round(ov$pct_en_replicados, 1), pct_resto = round(ov$pct_en_resto, 1), OR = round(ov$OR, 2), p = signif(ov$p_Fisher, 2), fdr = signif(ov$fdr, 2)))
} else S("  (omitido: faltan ", fp, " o ", fm, ")")

# ---- 5. Los cuatro genes de la version inicial (caso de estudio a posteriori) ---
SH("5. Genes de la version inicial del TFM, evaluados a posteriori con los mismos criterios")
g4 <- CANDIDATES[CANDIDATES %in% R_os$gen]
tab4 <- data.frame(gen = g4, z_SEQC = R_os$z_SEQC[match(g4, R_os$gen)], fdr_SEQC = R_os$fdr_SEQC[match(g4, R_os$gen)], p_SEQC = R_os$p_SEQC[match(g4, R_os$gen)],
                   z_TARGET = R_os$z_TARGET[match(g4, R_os$gen)], p_TARGET = R_os$p_TARGET[match(g4, R_os$gen)], replicado = R_os$replicado[match(g4, R_os$gen)])
PT(tab4); write.csv(tab4, OUT("tabla_4_genes_iniciales.csv"), row.names = FALSE)

# ---- 6. Figura -----------------------------------------------------------------
SH("6. Figura")
d <- R_os[is.finite(R_os$z_SEQC) & is.finite(R_os$z_TARGET), ]
d$cat <- ifelse(d$replicado & d$coef_SEQC > 0, "riesgo", ifelse(d$replicado, "prot", ifelse(d$descubierto, "solo", "resto")))
d$cat <- factor(d$cat, levels = c("resto", "solo", "prot", "riesgo"), labels = c("Resto", "Solo en SEQC", "Replicado, protecci\u00f3n", "Replicado, riesgo"))
lab <- d[d$gen %in% CANDIDATES, ]; lab$gen_f <- lab$gen
topl <- head(d[d$replicado, ][order(d$p_TARGET[d$replicado]), ], 8)
p <- ggplot(d[order(d$cat), ], aes(z_SEQC, z_TARGET, colour = cat)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(shape = cat, size = cat), alpha = 0.6) +
  scale_colour_manual(values = c("grey82", "grey35", "#2C7BB6", "#D7191C"), drop = FALSE, name = NULL) +
  scale_shape_manual(values = c(16, 16, 15, 17), drop = FALSE, name = NULL) + scale_size_manual(values = c(0.6, 0.9, 1.8, 1.8), drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.2)), size = "none") +
  geom_text_repel(data = lab, aes(label = gen), colour = "black", fontface = "italic", size = 3.4, max.overlaps = Inf, box.padding = 0.6) +
  geom_text_repel(data = topl, aes(label = gen), colour = "black", fontface = "italic", size = 3, max.overlaps = Inf) +
  labs(title = "Replicaci\u00f3n gen a gen del efecto sobre la supervivencia global",
       subtitle = sprintf("Estad\u00edstico z de Cox (ajustado por edad, amplificaci\u00f3n del oncog\u00e9n y estadio 4)\n%s genes comunes; Spearman entre cohortes = %s",
                          format(nrow(d), big.mark = ".", decimal.mark = ","), formatC(cor(d$z_SEQC, d$z_TARGET, method = "spearman"), format = "f", digits = 2, decimal.mark = ",")),
       x = sprintf("z en SEQC (n = %d)", nrow(sa)), y = sprintf("z en TARGET-NBL (n = %d)", nrow(st))) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot", plot.title = element_text(face = "bold", size = 13))
ggsave(OUT("fig_replica_genes.png"), p, width = 8, height = 6.2, dpi = 200, bg = "white"); S("  figura: fig_replica_genes.png")
end_log()
