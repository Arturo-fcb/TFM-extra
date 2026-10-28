# =============================================================================
# 03_GSE62564_cohorte_clinica.R   (cohorte SEQC, RNA-seq, 498 pacientes)
#   Parte A: expresion diferencial alto riesgo (HR) vs no alto riesgo, a nivel de GEN
#   Parte B: enriquecimiento funcional
#   Parte C: supervivencia (Kaplan-Meier, Cox) con comprobaciones de robustez
# Entradas : Datos/62564/GSE62564_series_matrix.txt  y  GSE62564_SEQC_NB_RNA-Seq_log2RPM(.txt/.txt.gz)
#            (opcional) Datos/49711/GSE49711_series_matrix.txt  -> estatus MYCN y estadio
# Salidas  : resultados/GSE62564/
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("limma", "survival", "AnnotationDbi", "org.Hs.eg.db", "ggplot2"))
suppressPackageStartupMessages({ library(limma); library(survival); library(AnnotationDbi); library(org.Hs.eg.db) })
start_log("GSE62564")
D <- file.path(DATA, "62564")
RUN_GENOMEWIDE <- TRUE      # cribado Cox de todo el transcriptoma (tarda unos minutos); poner FALSE para saltarlo

# -----------------------------------------------------------------------------
SH("1. Metadatos clinicos (series matrix)")
meta <- parse_series_matrix(file.path(D, "GSE62564_series_matrix.txt"))
meta$id <- sub(" \\[.*\\]$", "", meta$title)
cat("Columnas detectadas:", paste(names(meta), collapse = ", "), "\n")
num <- function(x) suppressWarnings(as.numeric(x))
meta$os_days   <- num(meta$os.day);  meta$os_event  <- num(meta$os.bin)
meta$efs_days  <- num(meta$efs.day); meta$efs_event <- num(meta$efs.bin)
meta$age_days  <- num(meta$age)
meta$HR        <- as.integer(!is.na(meta$high.risk) & meta$high.risk == "HR")   # solo los HR llevan la clave "high risk: HR"
meta$set       <- meta$training.validation
meta$sex       <- meta$Sex
rownames(meta) <- meta$id
S("Pacientes: ", nrow(meta), " | alto riesgo (HR): ", sum(meta$HR), " | no alto riesgo: ", sum(meta$HR == 0))
S("Eventos OS: ", sum(meta$os_event), " | eventos EFS: ", sum(meta$efs_event))
S("Training/validation: ", paste(names(table(meta$set)), table(meta$set), collapse = " / "))
S("Edad al diagnostico (dias): mediana ", median(meta$age_days), "; >=18 meses (547 d): ", sum(meta$age_days >= 547))

# -----------------------------------------------------------------------------
SH("2. Matriz de expresion y paso transcrito -> gen")
fx <- file.path(D, "GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt"); fz <- paste0(fx, ".gz")
cat("Leyendo matriz de expresion (puede tardar 1-2 min)...\n")
X <- if (file.exists(fx) && has("data.table")) { x <- data.table::fread(fx, data.table = FALSE); rownames(x) <- x[[1]]; as.matrix(x[, -1]) } else {
  src <- if (file.exists(fz)) gzfile(fz) else fx; as.matrix(read.delim(src, header = TRUE, row.names = 1, check.names = FALSE)) }
stopifnot(all(meta$id %in% colnames(X)))
X <- X[, meta$id]
S("Matriz: ", nrow(X), " filas (identificadores RefSeq = TRANSCRITOS) x ", ncol(X), " muestras")
S("Prefijos RefSeq: ", paste(names(table(substr(rownames(X), 1, 3))), table(substr(rownames(X), 1, 3)), collapse = "  "))
off <- 2^min(X)
S("Valor minimo de la matriz = ", round(min(X), 3), "  => pseudoconteo implicito RPM = ", signif(off, 3))

ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(X), keytype = "REFSEQ", columns = c("REFSEQ", "SYMBOL")))
ann <- ann[!is.na(ann$SYMBOL) & !duplicated(ann$REFSEQ), ]
S("Transcritos con simbolo: ", nrow(ann), " -> genes unicos: ", length(unique(ann$SYMBOL)))
iso_per_gene <- table(ann$SYMBOL); S("Genes con >1 transcrito en la matriz: ", sum(iso_per_gene > 1))
# expresion de gen = suma de RPM de sus transcritos (RPM es aditivo), de nuevo en log2
rpm <- pmax(2^X[ann$REFSEQ, ] - off, 0)
G <- log2(rowsum(rpm, ann$SYMBOL) + off); rm(rpm)
keepG <- rowSums(G >= 0) >= 0.2 * ncol(G)        # RPM>=1 en >=20% de las muestras
Gf <- G[keepG, ]
S("Genes a nivel gen: ", nrow(G), " | tras filtro de expresion (RPM>=1 en >=20% muestras): ", nrow(Gf))
cat("\nTranscritos RefSeq de los genes candidatos en la matriz (log2RPM medio):\n")
cand_iso <- ann[ann$SYMBOL %in% CANDIDATES, ]; cand_iso$mean_log2RPM <- rowMeans(X[cand_iso$REFSEQ, , drop = FALSE])
STAB(cand_iso[order(cand_iso$SYMBOL), ], 3)

# -----------------------------------------------------------------------------
SH("3. Expresion diferencial HR vs no-HR")
des <- model.matrix(~ HR, data = meta)
# 3a. Reproduccion del analisis ORIGINAL (a nivel de transcrito, sin filtro): para comparar con la memoria
tt_tx <- topTable(eBayes(lmFit(X, des)), coef = "HR", number = Inf)
dtx <- tt_tx[tt_tx$adj.P.Val < FDR_CUT & abs(tt_tx$logFC) > LFC_CUT, ]
S("[ORIGINAL, a nivel TRANSCRITO] 'DEGs' = ", nrow(dtx), " (up ", sum(dtx$logFC > 0), ", down ", sum(dtx$logFC < 0), ")   <- memoria: 3675 (1194 / 2481)")
S("   de esas filas, simbolos unicos: ", length(unique(ann$SYMBOL[ann$REFSEQ %in% rownames(dtx)])))
# 3b. Analisis CORREGIDO a nivel de GEN
fit <- eBayes(lmFit(Gf, des), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "HR", number = Inf)
write.csv(tt, OUT("DE_HR_vs_nonHR_gen.csv"))
S("[A NIVEL GEN] testados ", nrow(tt)); STAB(t(as.data.frame(deg_counts(tt))), 7)
S("Candidatos (logFC = HR vs no-HR; positivo = mas alto en HR):"); STAB(cand_table(tt), 3)
des2 <- model.matrix(~ HR + sex, data = meta)
tt_s <- topTable(eBayes(lmFit(Gf, des2), trend = TRUE, robust = TRUE), coef = "HR", number = Inf)
S("Sensibilidad ajustando por sexo: DEGs (FDR & |LFC|>1) = ", deg_counts(tt_s)["FDR_y_LFC"], " (sin ajustar: ", deg_counts(tt)["FDR_y_LFC"], ")")

volcano_plot(tt, "HR vs no-HR (GSE62564, nivel gen)", OUT("volcano_HR_vs_nonHR.png"), xlab = "log2FC (HR / no-HR)")
vv <- apply(Gf, 1, var); topv <- names(sort(vv, decreasing = TRUE))[1:5000]
pca_plot(Gf[topv, ], factor(ifelse(meta$HR == 1, "HR", "no-HR")), OUT("PCA_HR.png"), "PCA (5000 genes mas variables)")
topg <- rownames(tt)[1:50]
heatmap_png(Gf[topg, ], data.frame(Riesgo = factor(ifelse(meta$HR == 1, "HR", "no-HR")), row.names = meta$id), OUT("heatmap_top50.png"), "Top 50 genes (HR vs no-HR)")

SH("4. Enriquecimiento funcional (HR vs no-HR)")
run_enrichment(tt, "HRvsNonHR")

# -----------------------------------------------------------------------------
SH("5. Supervivencia: preparacion")
sv <- meta[, c("id", "os_days", "os_event", "efs_days", "efs_event", "age_days", "HR", "set", "sex")]
for (g in CANDIDATES) sv[[g]] <- as.numeric(G[g, sv$id])          # expresion a nivel gen (suma de transcritos)
for (g in CANDIDATES) sv[[paste0(g, "_sd")]] <- as.numeric(scale(sv[[g]]))
sv$agecat <- as.integer(sv$age_days >= 547)                         # >= 18 meses
# Transcritos usados en el analisis ORIGINAL (para reproducir la memoria)
TX_ORIG <- c(ARPP21 = "NM_016300", DLG2 = "NM_001364", AXL = "NM_001699", MYCN = "NM_005378")
for (g in CANDIDATES) sv[[paste0(g, "_tx")]] <- as.numeric(X[TX_ORIG[g], sv$id])
S("Correlacion gen-nivel vs transcrito-original (Spearman): ",
  paste(sprintf("%s %.3f", CANDIDATES, sapply(CANDIDATES, function(g) cor(sv[[g]], sv[[paste0(g, "_tx")]], method = "spearman"))), collapse = " | "))

cox1 <- function(form, d, term) {
  m <- tryCatch(coxph(form, data = d), error = function(e) NULL)
  if (is.null(m)) return(c(HR = NA, lo = NA, hi = NA, p = NA, n = NA, ev = NA, C = NA))
  s <- summary(m)
  c(HR = unname(s$conf.int[term, "exp(coef)"]), lo = unname(s$conf.int[term, "lower .95"]), hi = unname(s$conf.int[term, "upper .95"]),
    p = unname(s$coefficients[term, "Pr(>|z|)"]), n = s$n, ev = s$nevent, C = unname(s$concordance[1]))
}
coxtab <- function(rhs, d = sv, genes = CANDIDATES, time = "os_days", event = "os_event", var = function(g) g) {
  t(sapply(genes, function(g) cox1(as.formula(paste0("Surv(", time, ",", event, ") ~ ", var(g), rhs)), d, var(g))))
}

SH("6. Kaplan-Meier (corte en la mediana de cada gen)")
km_one <- function(d, g, time = "os_days", event = "os_event") {
  d$grp <- factor(ifelse(d[[g]] >= median(d[[g]]), "Alto", "Bajo"), levels = c("Bajo", "Alto"))
  f <- as.formula(paste0("Surv(", time, "/365.25,", event, ") ~ grp"))
  fit <- survfit(f, data = d); lr <- survdiff(f, data = d)
  s5 <- summary(fit, times = 5, extend = TRUE)$surv; s10 <- summary(fit, times = 10, extend = TRUE)$surv
  c(n_bajo = sum(d$grp == "Bajo"), n_alto = sum(d$grp == "Alto"), S5_bajo = s5[1], S5_alto = s5[2], S10_bajo = s10[1], S10_alto = s10[2],
    logrank_p = 1 - pchisq(lr$chisq, 1), pct_HR_bajo = 100 * mean(d$HR[d$grp == "Bajo"]), pct_HR_alto = 100 * mean(d$HR[d$grp == "Alto"]))
}
S("OS, toda la cohorte (n=498):");           km_all <- t(sapply(CANDIDATES, function(g) km_one(sv, g)));            STAB(km_all, 3)
S("EFS, toda la cohorte:");                  STAB(t(sapply(CANDIDATES, function(g) km_one(sv, g, "efs_days", "efs_event"))), 3)
S("OS solo no-HR (mediana recalculada dentro del grupo):"); STAB(t(sapply(CANDIDATES, function(g) km_one(sv[sv$HR == 0, ], g))), 3)
S("OS solo HR (mediana recalculada dentro del grupo):");    STAB(t(sapply(CANDIDATES, function(g) km_one(sv[sv$HR == 1, ], g))), 3)
S("Memoria: ARPP21 ~92% vs ~60%, DLG2 88/65, AXL 83/67, MYCN 68/82 (supervivencia a 10 anios, alto/bajo)")
write.csv(km_all, OUT("KM_resumen_OS.csv"))
for (g in CANDIDATES) {
  d <- sv; d$grp <- factor(ifelse(d[[g]] >= median(d[[g]]), paste(g, "alto"), paste(g, "bajo")))
  d$t <- d$os_days / 365.25
  fpng <- OUT(paste0("KM_OS_", g, ".png")); ok <- FALSE
  if (has("survminer")) ok <- tryCatch({                      # grafico bonito; si survminer falla (p. ej. conflicto de versiones), se usa el basico
    fit <- survfit(Surv(t, os_event) ~ grp, data = d)
    p <- survminer::ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE, xlab = "Tiempo (anios)", ylab = "Supervivencia global", title = g)
    png(fpng, 1800, 1500, res = 200); print(p); dev.off(); TRUE
  }, error = function(e) { while (dev.cur() > 1) dev.off(); message("  (survminer fallo para ", g, ": ", substr(conditionMessage(e), 1, 90), " -> grafico basico)"); FALSE })
  if (!ok) {
    png(fpng, 1500, 1200, res = 200)
    plot(survfit(Surv(t, os_event) ~ grp, data = d), col = c("red", "blue"), lwd = 2, xlab = "Tiempo (anios)", ylab = "Supervivencia global", main = g)
    legend("bottomleft", legend = levels(d$grp), col = c("red", "blue"), lwd = 2); dev.off()
  }
}

SH("7. Regresion de Cox")
S("(HR por unidad de log2 RPM; 'C' = indice de concordancia)")
S("7a. Univariante, OS");                             u <- coxtab("");                           STAB(u, 3)
S("7b. Univariante, HR por 1 desviacion estandar");   STAB(coxtab("", var = function(g) paste0(g, "_sd")), 3)
S("7c. ORIGINAL de la memoria: gen + edad (dias) + alto riesgo, OS (n=498, 105 eventos)")
m2 <- coxtab(" + age_days + HR"); m2 <- cbind(m2, p_BH4 = p.adjust(m2[, "p"], "BH")); STAB(m2, 4)
S("   esperado en memoria: ARPP21 0.882 (0.803-0.970) p=0.009 | DLG2 0.940 p=0.317 | AXL 0.870 p=0.073 | MYCN 1.162 (1.063-1.270) p=0.001")
S("7d. Igual pero con el TRANSCRITO original (NM_016300, NM_001364, NM_001699, NM_005378) en vez del gen: debe reproducir 7c")
STAB(coxtab(" + age_days + HR", var = function(g) paste0(g, "_tx")), 4)
S("7e. Edad como ESTRATO (<18 / >=18 meses) en vez de covariable (la edad incumple riesgos proporcionales)")
STAB(coxtab(" + strata(agecat) + HR"), 4)
S("7f. Edad en logaritmo + HR");                      STAB(coxtab(" + log(age_days + 30) + HR"), 4)
S("7g. Solo mitad TRAINING de SEQC (gen + edad + HR)");   STAB(coxtab(" + age_days + HR", d = sv[sv$set == "Training", ]), 4)
S("7h. Solo mitad VALIDATION de SEQC (gen + edad + HR)"); STAB(coxtab(" + age_days + HR", d = sv[sv$set == "Validation", ]), 4)
S("7i. Dentro de no-HR (gen + edad); OJO pocos eventos: ", sum(sv$os_event[sv$HR == 0]));  STAB(coxtab(" + age_days", d = sv[sv$HR == 0, ]), 4)
S("7j. Dentro de HR (gen + edad), eventos: ", sum(sv$os_event[sv$HR == 1]));               STAB(coxtab(" + age_days", d = sv[sv$HR == 1, ]), 4)
S("7k. Interaccion gen x HR: p del test de razon de verosimilitud")
inter <- sapply(CANDIDATES, function(g) { a <- coxph(as.formula(paste0("Surv(os_days,os_event) ~ ", g, " + age_days + HR")), data = sv)
  b <- coxph(as.formula(paste0("Surv(os_days,os_event) ~ ", g, " * HR + age_days")), data = sv); anova(a, b)$`Pr(>|Chi|)`[2] }); STAB(t(inter), 4)
S("7l. Mismo modelo 7c para EFS (supervivencia libre de evento)")
STAB(coxtab(" + age_days + HR", time = "efs_days", event = "efs_event"), 4)

SH("8. Riesgos proporcionales (cox.zph) del modelo 7c")
zph <- lapply(CANDIDATES, function(g) cox.zph(coxph(as.formula(paste0("Surv(os_days,os_event) ~ ", g, " + age_days + HR")), data = sv)))
names(zph) <- CANDIDATES
for (g in CANDIDATES) { z <- zph[[g]]$table; S(sprintf("  %-7s p(gen)=%.3g p(edad)=%.3g p(HR)=%.3g p(GLOBAL)=%.3g", g, z[g, "p"], z["age_days", "p"], z["HR", "p"], z["GLOBAL", "p"])) }

SH("9. Correlaciones (Spearman) entre genes, edad y alto riesgo")
STAB(round(cor(sv[, c(CANDIDATES, "age_days", "HR")], method = "spearman"), 2), 2)

SH("10. Modelo ampliado con estatus MYCN y estadio (GSE49711, mismos pacientes)")
f49 <- file.path(DATA, "49711", "GSE49711_series_matrix.txt")
if (file.exists(f49)) {
  m49 <- parse_series_matrix(f49); m49$id <- sub(" \\[.*\\]$", "", m49$title)
  S("Claves clinicas en GSE49711: ", paste(setdiff(names(m49), c("title", "geo_accession", "source", "id")), collapse = ", "))
  kmy <- grep("mycn", names(m49), ignore.case = TRUE, value = TRUE); kst <- grep("stage", names(m49), ignore.case = TRUE, value = TRUE)
  S("Columnas candidatas MYCN: ", paste(kmy, collapse = ","), " | estadio: ", paste(kst, collapse = ","))
  ov <- intersect(sv$id, m49$id); S("Muestras en comun con GSE62564: ", length(ov))
  if (length(ov) > 400 && length(kmy) && length(kst)) {
    m49 <- m49[match(sv$id, m49$id), ]
    sv$mycn_amp <- as.integer(grepl("amplif", m49[[kmy[1]]], ignore.case = TRUE) & !grepl("non|not", m49[[kmy[1]]], ignore.case = TRUE))
    sv$stage4   <- as.integer(grepl("4", m49[[kst[1]]]) & !grepl("4s", m49[[kst[1]]], ignore.case = TRUE))
    print(table(MYCN = m49[[kmy[1]]], useNA = "ifany")); print(table(Estadio = m49[[kst[1]]], useNA = "ifany"))
    S("MYCN amplificado: ", sum(sv$mycn_amp, na.rm = TRUE), " | estadio 4: ", sum(sv$stage4, na.rm = TRUE))
    S("10a. gen + edad + MYCN_amp + estadio4 (SIN la clasificacion de riesgo, que ya incluye esas variables)")
    STAB(coxtab(" + age_days + mycn_amp + stage4"), 4)
    S("10b. gen + edad + MYCN_amp + estadio4 + HR");          STAB(coxtab(" + age_days + mycn_amp + stage4 + HR"), 4)
    S("10c. Dentro de MYCN no amplificado (gen + edad + estadio4)")
    STAB(coxtab(" + age_days + stage4", d = sv[which(sv$mycn_amp == 0), ]), 4)
  } else S("  No se pudo construir MYCN/estadio automaticamente: revisa las claves impresas arriba y avisame.")
} else S("  (omitido) Falta Datos/49711/GSE49711_series_matrix.txt. Descargalo de GEO (GSE49711 > 'Series Matrix File(s)') y vuelve a ejecutar.")

SH("11. Cribado de TODO el transcriptoma con el mismo modelo (gen + edad + HR): contexto para ARPP21")
if (RUN_GENOMEWIDE) {
  cat("Ajustando ", nrow(Gf), " modelos de Cox (varios minutos)...\n")
  res <- t(vapply(seq_len(nrow(Gf)), function(i) {
    v <- Gf[i, ]; m <- tryCatch(summary(coxph(Surv(sv$os_days, sv$os_event) ~ v + sv$age_days + sv$HR))$coefficients[1, ], error = function(e) rep(NA_real_, 5))
    c(coef = m[1], p = m[5]) }, numeric(2)))
  rownames(res) <- rownames(Gf); res <- as.data.frame(res); colnames(res) <- c("coef", "p"); res$fdr <- p.adjust(res$p, "BH")
  res$rank_p <- rank(res$p); write.csv(res[order(res$p), ], OUT("Cox_transcriptoma_completo.csv"))
  S("Genes testados: ", nrow(res), " | con p<0.05: ", sum(res$p < 0.05, na.rm = TRUE), " | con FDR<0.05: ", sum(res$fdr < 0.05, na.rm = TRUE))
  S("Posicion de cada candidato por p-valor:"); STAB(res[CANDIDATES, c("coef", "p", "fdr", "rank_p")], 3)
} else S("  (saltado: RUN_GENOMEWIDE = FALSE)")

end_log()
