# =============================================================================
# 04_GSE62564_ajuste_MYCN_estadio.R
#   Misma cohorte SEQC (GSE62564, 498 pacientes) + variables clinicas de GSE49711
#   (estatus MYCN, estadio INSS) para los MISMOS pacientes (SEQC_NBxxx).
#   Analisis SIN genes candidatos a priori: todo el transcriptoma.
#     1. Pacientes comunes y coherencia entre GSE62564 y GSE49711
#     2. Cuanto "explican" edad/MYCN/estadio la etiqueta de alto riesgo (HR)
#     3. Que parte de la firma HR vs no-HR sobrevive al ajustar por estatus MYCN
#     4. Cribado Cox de todo el transcriptoma con edad + MYCN + estadio (sin HR) y con HR
#     5. Enriquecimiento funcional del cribado Cox ajustado
#     6. Genes concordantes entre datasets: ?sobreviven al ajuste?
# Entradas : Datos/62564/*  y  Datos/49711/GSE49711_series_matrix.txt(.gz)  (o Datos/GSE49711_series_matrix.txt.gz)
# Salidas  : resultados/GSE62564_ajuste/
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("limma", "survival", "AnnotationDbi", "org.Hs.eg.db"))
suppressPackageStartupMessages({ library(limma); library(survival); library(AnnotationDbi); library(org.Hs.eg.db) })
start_log("GSE62564_ajuste")
D <- file.path(DATA, "62564")
num <- function(x) suppressWarnings(as.numeric(x))
PT <- function(x) for (l in capture.output(print(x))) S(l)      # imprime cualquier objeto tambien en el RESUMEN

# -----------------------------------------------------------------------------
SH("1. Metadatos GSE62564 + matriz de expresion a nivel gen (igual que el script 03)")
meta <- parse_series_matrix(file.path(D, "GSE62564_series_matrix.txt"))
meta$id <- sub(" \\[.*\\]$", "", meta$title)
meta$os_days <- num(meta$os.day); meta$os_event <- num(meta$os.bin)
meta$efs_days <- num(meta$efs.day); meta$efs_event <- num(meta$efs.bin)
meta$age_days <- num(meta$age)
meta$HR <- as.integer(!is.na(meta$high.risk) & meta$high.risk == "HR")
meta$set <- meta$training.validation
rownames(meta) <- meta$id

cache <- file.path(RES, "cache_G_GSE62564.rds")
if (file.exists(cache)) { G <- readRDS(cache); S("Matriz gen leida de cache: ", nrow(G), " genes") } else {
  fx <- file.path(D, "GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt"); fz <- paste0(fx, ".gz")
  cat("Leyendo matriz de expresion (1-2 min)...\n")
  X <- if (file.exists(fx) && has("data.table")) { x <- data.table::fread(fx, data.table = FALSE); rownames(x) <- x[[1]]; as.matrix(x[, -1]) } else {
    src <- if (file.exists(fz)) gzfile(fz) else fx; as.matrix(read.delim(src, header = TRUE, row.names = 1, check.names = FALSE)) }
  X <- X[, meta$id]; off <- 2^min(X)
  ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(X), keytype = "REFSEQ", columns = c("REFSEQ", "SYMBOL")))
  ann <- ann[!is.na(ann$SYMBOL) & !duplicated(ann$REFSEQ), ]
  rpm <- pmax(2^X[ann$REFSEQ, ] - off, 0); G <- log2(rowsum(rpm, ann$SYMBOL) + off); rm(rpm, X)
  saveRDS(G, cache); S("Matriz gen calculada: ", nrow(G), " genes")
}
Gf <- G[rowSums(G >= 0) >= 0.2 * ncol(G), meta$id]
S("Genes tras filtro (RPM>=1 en >=20% muestras): ", nrow(Gf), " | pacientes: ", ncol(Gf))

# -----------------------------------------------------------------------------
SH("2. Variables clinicas de GSE49711 (mismos pacientes) y coherencia con GSE62564")
cand <- c(file.path(DATA, "49711", "GSE49711_series_matrix.txt"), file.path(DATA, "49711", "GSE49711_series_matrix.txt.gz"),
          file.path(DATA, "GSE49711_series_matrix.txt"), file.path(DATA, "GSE49711_series_matrix.txt.gz"))
f49 <- cand[file.exists(cand)][1]
if (is.na(f49)) stop("No encuentro GSE49711_series_matrix.txt(.gz) en Datos/ ni en Datos/49711/")
S("Archivo GSE49711: ", f49)
m49 <- parse_series_matrix(f49); m49$id <- sub(" \\[.*\\]$", "", m49$title)
S("Claves clinicas GSE49711: ", paste(setdiff(names(m49), c("title", "geo_accession", "source", "id")), collapse = ", "))
S("Pacientes comunes: ", length(intersect(meta$id, m49$id)), " de ", nrow(meta))
m49 <- m49[match(meta$id, m49$id), ]
sv <- meta[, c("id", "os_days", "os_event", "efs_days", "efs_event", "age_days", "HR", "set")]
sv$mycn_amp <- ifelse(m49$mycn.status %in% c("0", "1"), as.integer(m49$mycn.status), NA)
sv$stage    <- factor(sub("^\\s+|\\s+$", "", m49$inss.stage), levels = c("1", "2", "3", "4S", "4"))
sv$stage4   <- as.integer(sv$stage == "4")
sv$HR_49711 <- ifelse(m49$high.risk %in% c("0", "1"), as.integer(m49$high.risk), NA)
sv$age49711 <- num(m49$age.at.diagnosis)
S("Coherencia HR (GSE62564 vs GSE49711):"); PT(table(HR_62564 = sv$HR, HR_49711 = sv$HR_49711, useNA = "ifany"))
S("Coherencia edad: correlacion = ", round(cor(sv$age_days, sv$age49711, use = "complete.obs"), 4), " | iguales en ", sum(sv$age_days == sv$age49711, na.rm = TRUE), " de ", nrow(sv))
S("Estatus MYCN: ");  PT(table(MYCN = sv$mycn_amp, useNA = "ifany"))
S("Estadio INSS: ");  PT(table(sv$stage, useNA = "ifany"))
sv <- sv[!is.na(sv$mycn_amp) & !is.na(sv$stage), ]
rownames(sv) <- sv$id
S("Pacientes con MYCN y estadio conocidos: ", nrow(sv), " | MYCN amplificado: ", sum(sv$mycn_amp), " | estadio 4: ", sum(sv$stage4))
sv$logage <- log1p(sv$age_days)
sv$age18  <- as.integer(sv$age_days >= 547)

# -----------------------------------------------------------------------------
SH("3. ?Cuanto determinan edad, MYCN y estadio la etiqueta de alto riesgo (HR)?")
S("HR segun MYCN y estadio (y edad >=18 meses):")
PT(ftable(MYCN = sv$mycn_amp, Estadio = sv$stage, Edad18m = sv$age18, HR = sv$HR))
S("Tabla HR por estatus MYCN:"); PT(table(MYCN = sv$mycn_amp, HR = sv$HR))
glm_hr <- glm(HR ~ mycn_amp + stage4 + age18, data = sv, family = binomial)
S("Prediccion de HR con MYCN + estadio4 + edad>=18m: exactitud = ",
  round(mean((fitted(glm_hr) > 0.5) == (sv$HR == 1)), 3))
S("Eventos OS: ", sum(sv$os_event), " | eventos EFS: ", sum(sv$efs_event))

# -----------------------------------------------------------------------------
SH("4. Firma HR vs no-HR: ?que parte sobrevive al ajustar por estatus MYCN?")
Gs <- Gf[, sv$id]
lim <- function(form, d = sv, M = Gs, coef) {
  des <- model.matrix(form, data = d)
  tt <- topTable(eBayes(lmFit(M, des), trend = TRUE, robust = TRUE), coef = coef, number = Inf)
  tt
}
tt_hr0  <- lim(~ HR, coef = "HR")
tt_hr_m <- lim(~ HR + mycn_amp, coef = "HR")
tt_mycn <- lim(~ mycn_amp, coef = "mycn_amp")
ix <- sv$mycn_amp == 0
tt_hr_nonamp <- lim(~ HR, d = sv[ix, ], M = Gs[, ix], coef = "HR")
tt_stage4 <- lim(~ stage4 + mycn_amp, coef = "stage4")
cnt <- rbind(`HR vs no-HR (sin ajustar)` = deg_counts(tt_hr0), `HR vs no-HR ajustado por MYCN` = deg_counts(tt_hr_m),
             `MYCN amp vs no amp (todos)` = deg_counts(tt_mycn), `HR vs no-HR solo en MYCN no amp` = deg_counts(tt_hr_nonamp),
             `estadio 4 vs resto, ajustado por MYCN` = deg_counts(tt_stage4))
STAB(cnt, 7)
cm <- intersect(rownames(tt_hr0), rownames(tt_mycn))
S("Correlacion Spearman de t: HR sin ajustar vs MYCN amp = ", round(cor(tt_hr0[cm, "t"], tt_mycn[cm, "t"], method = "spearman"), 3),
  " | HR sin ajustar vs HR ajustado por MYCN = ", round(cor(tt_hr0[cm, "t"], tt_hr_m[cm, "t"], method = "spearman"), 3),
  " | HR sin ajustar vs HR solo MYCN no amp = ", round(cor(tt_hr0[cm, "t"], tt_hr_nonamp[cm, "t"], method = "spearman"), 3))
s0 <- rownames(tt_hr0)[tt_hr0$adj.P.Val < FDR_CUT & abs(tt_hr0$logFC) > LFC_CUT]
s1 <- rownames(tt_hr_m)[tt_hr_m$adj.P.Val < FDR_CUT & abs(tt_hr_m$logFC) > LFC_CUT]
S("DEGs HR (FDR<0.05, |log2FC|>1): ", length(s0), " sin ajustar | ", length(s1), " ajustando por MYCN | comunes: ", length(intersect(s0, s1)))
write.csv(tt_hr0, OUT("DE_HR_sin_ajustar_n_comun.csv")); write.csv(tt_hr_m, OUT("DE_HR_ajustado_MYCN.csv"))
write.csv(tt_mycn, OUT("DE_MYCNamp_vs_noamp.csv")); write.csv(tt_hr_nonamp, OUT("DE_HR_en_MYCN_noamp.csv"))
S("Top 15 genes por FDR, MYCN amp vs no amp (logFC>0 = mas alto en amplificados):")
STAB(round(head(tt_mycn[, c("logFC", "AveExpr", "t", "adj.P.Val")], 15), 4), 4)
S("Top 15 genes por FDR, HR vs no-HR ajustado por MYCN:")
STAB(round(head(tt_hr_m[, c("logFC", "AveExpr", "t", "adj.P.Val")], 15), 4), 4)

if (has("clusterProfiler")) { S("-- Enriquecimiento HR vs no-HR AJUSTADO por MYCN --"); run_enrichment(tt_hr_m, "HR_ajustado_MYCN")
                              S("-- Enriquecimiento MYCN amp vs no amp --"); run_enrichment(tt_mycn, "MYCNamp_vs_noamp") }

# -----------------------------------------------------------------------------
SH("5. Cribado Cox de TODO el transcriptoma (OS), con distintos ajustes")
genomewide <- function(rhs, d, M, time = "os_days", event = "os_event") {
  M <- M[, d$id]; dd <- d
  res <- t(vapply(seq_len(nrow(M)), function(i) {
    dd$v <- M[i, ]
    m <- tryCatch(suppressWarnings(summary(coxph(as.formula(paste0("Surv(", time, ",", event, ") ~ v", rhs)), data = dd))$coefficients[1, ]),
                  error = function(e) rep(NA_real_, 5))
    c(m[1], m[5]) }, numeric(2)))
  res <- as.data.frame(res); colnames(res) <- c("coef", "p"); rownames(res) <- rownames(M)
  res$z <- sign(res$coef) * qnorm(res$p / 2, lower.tail = FALSE); res$fdr <- p.adjust(res$p, "BH"); res$rank_p <- rank(res$p)
  res
}
cat("Ajustando modelos de Cox genoma completo (varios minutos por modelo)...\n")
t0 <- Sys.time()
cx <- list()
cx$M0_edad_HR        <- genomewide(" + logage + HR", sv, Gs)                       # referencia: modelo del script 03 (con log edad)
cx$M1_edad_MYCN_est  <- genomewide(" + logage + mycn_amp + stage4", sv, Gs)       # SIN HR
cx$M2_edad_MYCN_est_HR <- genomewide(" + logage + mycn_amp + stage4 + HR", sv, Gs)
cx$M3_MYCNnoamp      <- genomewide(" + logage + stage4", sv[sv$mycn_amp == 0, ], Gs)
cx$M4_EFS_edad_MYCN_est <- genomewide(" + logage + mycn_amp + stage4", sv, Gs, "efs_days", "efs_event")
S("Tiempo de calculo: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
for (n in names(cx)) write.csv(cx[[n]][order(cx[[n]]$p), ], OUT(paste0("Cox_", n, ".csv")))
ev <- c(M0_edad_HR = sum(sv$os_event), M1_edad_MYCN_est = sum(sv$os_event), M2_edad_MYCN_est_HR = sum(sv$os_event),
        M3_MYCNnoamp = sum(sv$os_event[sv$mycn_amp == 0]), M4_EFS_edad_MYCN_est = sum(sv$efs_event))
nn <- c(M0_edad_HR = nrow(sv), M1_edad_MYCN_est = nrow(sv), M2_edad_MYCN_est_HR = nrow(sv), M3_MYCNnoamp = sum(sv$mycn_amp == 0), M4_EFS_edad_MYCN_est = nrow(sv))
sm <- t(sapply(names(cx), function(n) c(pacientes = nn[[n]], eventos = ev[[n]], genes = nrow(cx[[n]]), p_menor_0.05 = sum(cx[[n]]$p < 0.05, na.rm = TRUE),
                                        FDR_menor_0.05 = sum(cx[[n]]$fdr < 0.05, na.rm = TRUE), FDR_menor_0.05_riesgo = sum(cx[[n]]$fdr < 0.05 & cx[[n]]$coef > 0, na.rm = TRUE),
                                        FDR_menor_0.05_proteccion = sum(cx[[n]]$fdr < 0.05 & cx[[n]]$coef < 0, na.rm = TRUE))))
S("Resumen del cribado (coef = log HR por unidad de log2 RPM; riesgo = coef>0):"); STAB(sm, 7)
S("Correlacion Spearman de z entre modelos:")
cmn <- Reduce(intersect, lapply(cx, rownames)); Z <- sapply(cx, function(r) r[cmn, "z"]); STAB(round(cor(Z, method = "spearman", use = "pairwise"), 2), 2)
S("Top 20 genes por p en M1 (edad + MYCN + estadio, sin HR):")
STAB(round(head(cx$M1_edad_MYCN_est[order(cx$M1_edad_MYCN_est$p), c("coef", "p", "fdr")], 20), 6), 5)
S("Top 20 genes por p en M2 (edad + MYCN + estadio + HR):")
STAB(round(head(cx$M2_edad_MYCN_est_HR[order(cx$M2_edad_MYCN_est_HR$p), c("coef", "p", "fdr")], 20), 6), 5)
S("Top 20 genes por p en M3 (solo MYCN no amplificado):")
STAB(round(head(cx$M3_MYCNnoamp[order(cx$M3_MYCNnoamp$p), c("coef", "p", "fdr")], 20), 6), 5)

SH("6. Enriquecimiento funcional del cribado Cox ajustado (ranking por z)")
for (n in c("M1_edad_MYCN_est", "M2_edad_MYCN_est_HR")) {
  r <- cx[[n]]; r <- r[is.finite(r$z), ]
  tt <- data.frame(logFC = r$coef, t = r$z, adj.P.Val = r$fdr, row.names = rownames(r))
  S("-- ", n, " (ORA: genes de riesgo/proteccion con FDR<0.05; GSEA por z) --")
  run_enrichment(tt, n, lfc = 0, kegg = FALSE)
}

# -----------------------------------------------------------------------------
SH("7. Genes concordantes entre datasets (definidos SIN mirar supervivencia ajustada nueva): ?sobreviven al ajuste?")
f_g <- file.path(RES, "GSE18271", "DE_GNB_vs_GN_gen.csv"); f_h <- file.path(RES, "GSE62564", "DE_HR_vs_nonHR_gen.csv")
f_x <- file.path(RES, "GSE62564", "Cox_transcriptoma_completo.csv"); f_c <- file.path(RES, "GSE212877", "D_pareado_mixto.csv")
if (all(file.exists(c(f_g, f_h, f_x, f_c)))) {
  g <- read.csv(f_g, row.names = 1); h <- read.csv(f_h, row.names = 1); x <- read.csv(f_x, row.names = 1); cl <- read.csv(f_c, row.names = 1)
  co <- Reduce(intersect, list(rownames(g), rownames(h), rownames(x), rownames(cl)))
  riesgo <- co[g[co, "adj.P.Val"] < .05 & h[co, "adj.P.Val"] < .05 & x[co, "fdr"] < .05 & g[co, "logFC"] > 0 & h[co, "logFC"] > 0 & x[co, "coef"] > 0]
  prot   <- co[g[co, "adj.P.Val"] < .05 & h[co, "adj.P.Val"] < .05 & x[co, "fdr"] < .05 & g[co, "logFC"] < 0 & h[co, "logFC"] < 0 & x[co, "coef"] < 0]
  S("Conjunto 'riesgo' (GNB>GN, HR>noHR, Cox gen+edad+HR coef>0; los tres FDR<0.05): ", length(riesgo), " genes | 'proteccion': ", length(prot))
  chk <- function(set, sign) {
    out <- sapply(c("M0_edad_HR", "M1_edad_MYCN_est", "M2_edad_MYCN_est_HR", "M3_MYCNnoamp", "M4_EFS_edad_MYCN_est"), function(n) {
      r <- cx[[n]]; s <- intersect(set, rownames(r)); c(en_modelo = length(s), FDR05_mismo_signo = sum(r[s, "fdr"] < .05 & sign * r[s, "coef"] > 0, na.rm = TRUE),
                                                        p05_mismo_signo = sum(r[s, "p"] < .05 & sign * r[s, "coef"] > 0, na.rm = TRUE), mismo_signo = sum(sign * r[s, "coef"] > 0, na.rm = TRUE)) })
    t(out) }
  S("Conjunto 'riesgo': cuantos mantienen el efecto en cada modelo"); STAB(chk(riesgo, +1), 6)
  S("Conjunto 'proteccion': cuantos mantienen el efecto en cada modelo"); STAB(chk(prot, -1), 6)
  write.csv(data.frame(gen = c(riesgo, prot), conjunto = c(rep("riesgo", length(riesgo)), rep("proteccion", length(prot)))), OUT("genes_concordantes_3_datasets.csv"), row.names = FALSE)
  # ?estan los genes concordantes sobre-representados entre los DEG de MYCN?
  mm <- intersect(rownames(tt_mycn), co)
  up_m <- mm[tt_mycn[mm, "adj.P.Val"] < .05 & tt_mycn[mm, "logFC"] > 0]
  tabl <- matrix(c(sum(riesgo %in% up_m), sum(!riesgo %in% up_m), sum(setdiff(mm, riesgo) %in% up_m), sum(!setdiff(mm, riesgo) %in% up_m)), 2)
  resto <- setdiff(mm, riesgo)
  S("Genes 'riesgo' que ademas suben en MYCN amplificado (FDR<0.05): ", tabl[1, 1], " de ", length(riesgo), " (", round(100 * tabl[1, 1] / length(riesgo)), "%) vs ",
    round(100 * sum(resto %in% up_m) / length(resto)), "% en el resto de genes | Fisher p = ", signif(fisher.test(tabl)$p.value, 3))
} else S("  (omitido: faltan los CSV de resultados de los scripts 01-03)")

end_log()
