# =============================================================================
# 06_validacion_TARGET_NBL.R   VALIDACION EXTERNA (cohorte independiente: TARGET-NBL, RNA-seq STAR del GDC)
#   Se aplican las MISMAS definiciones de programas (genes_programas.csv, creado por el script 05)
#   y el MISMO modelo Cox que en SEQC. No se reajusta nada en la cohorte nueva.
# Entradas : Datos/GDC NBL/  (o Datos/TARGET_NBL/)  con
#              - el "Sample Sheet" del GDC (TSV; nombre que contenga 'sample' y termine en .tsv o .txt)
#              - clinical.cohort*.json   (clinico del GDC, TARGET-NBL)
#              - los archivos *.rna_seq.augmented_star_gene_counts.tsv (en esa carpeta o en subcarpetas)
#            resultados/GSE62564_score/genes_programas.csv  (script 05)
# Salidas  : resultados/TARGET_NBL/
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("survival", "jsonlite"))
suppressPackageStartupMessages({ library(survival); library(jsonlite) })
start_log("TARGET_NBL")
# Carpetas donde buscar (todo dentro de extraordinaria/Datos): Datos/GDC NBL (sample sheet), Datos/GDC NBL/GDC (JSON clinico y carpetas de conteos)
# o, alternativamente, Datos/TARGET_NBL
TD <- file.path(DATA, "GDC NBL")
SEARCH <- unique(c(TD, file.path(TD, "GDC"), file.path(DATA, "TARGET_NBL"), DATA)); SEARCH <- SEARCH[dir.exists(SEARCH)]
COUNT_DIRS <- unique(c(TD, file.path(DATA, "TARGET_NBL"))); COUNT_DIRS <- COUNT_DIRS[dir.exists(COUNT_DIRS)]
S("Carpetas de busqueda: ", paste(SEARCH, collapse = " | "))
find1 <- function(pat) { for (d in SEARCH) { f <- list.files(d, full.names = TRUE); f <- f[grepl(pat, basename(f), ignore.case = TRUE)]; if (length(f)) return(f[1]) }; NA_character_ }
num <- function(x) suppressWarnings(as.numeric(x))
PT <- function(x) for (l in capture.output(print(x))) S(l)
nz <- function(x) if (is.null(x) || !length(x)) NA else x

# -----------------------------------------------------------------------------
SH("1. Clinica TARGET-NBL (JSON del GDC)")
fj <- find1("^clinical.*\\.json$")
if (is.na(fj)) stop("No encuentro clinical*.json en: ", paste(SEARCH, collapse = ", "))
S("JSON clinico: ", fj)
js <- fromJSON(fj, simplifyVector = FALSE)
S("Casos en el JSON: ", length(js))
cl <- do.call(rbind, lapply(js, function(c) {
  dg <- if (length(c$diagnoses)) c$diagnoses[[1]] else list()
  dm <- c$demographic; fu <- c$follow_ups
  myc <- unlist(lapply(fu, function(f) vapply(f$molecular_tests, function(m) if (identical(m$gene_symbol, "MYCN")) m$test_result else NA_character_, "")))
  myc <- myc[!is.na(myc)]
  mycn <- if ("Amplified" %in% myc) 1 else if ("Not Amplified" %in% myc) 0 else NA
  ev <- unlist(lapply(fu, function(f) if (!is.null(f$first_event) && !f$first_event %in% c("Censored")) f$days_to_first_event else NULL))
  evn <- unlist(lapply(fu, function(f) if (!is.null(f$first_event) && !f$first_event %in% c("Censored")) f$first_event else NULL))
  data.frame(case = c$submitter_id, age_days = num(nz(dg$age_at_diagnosis)), stage = nz(dg$inss_stage),
             risk_cog = nz(dg$cog_neuroblastoma_risk_group), vital = nz(dm$vital_status), days_death = num(nz(dm$days_to_death)),
             days_lfu = num(nz(dg$days_to_last_follow_up)), mycn_amp = mycn,
             efs_ev_days = if (length(ev)) min(num(ev), na.rm = TRUE) else NA_real_, n_ev = length(evn), stringsAsFactors = FALSE)
}))
cl <- cl[!is.na(cl$case) & !duplicated(cl$case), ]
cl$os_event <- as.integer(cl$vital == "Dead")
cl$os_days  <- ifelse(cl$os_event == 1 & !is.na(cl$days_death), cl$days_death, cl$days_lfu)
cl$efs_event <- as.integer(!is.na(cl$efs_ev_days) | cl$os_event == 1)
cl$efs_days <- ifelse(!is.na(cl$efs_ev_days), pmin(cl$efs_ev_days, ifelse(is.na(cl$os_days), Inf, cl$os_days)), cl$os_days)
cl$stage4 <- as.integer(cl$stage == "Stage 4")
cl$HR_cog <- as.integer(cl$risk_cog == "High Risk")
S("Casos con clinica: ", nrow(cl), " | MYCN conocido: ", sum(!is.na(cl$mycn_amp)), " (amplificado ", sum(cl$mycn_amp == 1, na.rm = TRUE), ")")
S("Estadio INSS:"); PT(table(cl$stage, useNA = "ifany")); S("Riesgo COG:"); PT(table(cl$risk_cog, useNA = "ifany"))

# -----------------------------------------------------------------------------
SH("2. Muestras de RNA-seq (Sample Sheet) y matriz de expresion")
ss_f <- find1("sample.*\\.(tsv|txt)$")
if (is.na(ss_f)) stop("No encuentro el Sample Sheet (TSV) en: ", paste(SEARCH, collapse = ", "))
ss <- read.delim(ss_f, check.names = FALSE, stringsAsFactors = FALSE)
S("Sample sheet: ", basename(ss_f), " | filas: ", nrow(ss), " | columnas: ", paste(names(ss), collapse = ", "))
colcase <- grep("^Case ID", names(ss), value = TRUE)[1]; colfile <- grep("^File Name", names(ss), value = TRUE)[1]
colsty <- grep("^Sample Type", names(ss), value = TRUE)[1]; if (is.na(colsty)) colsty <- grep("^Tumor Descriptor", names(ss), value = TRUE)[1]
colsid <- grep("^Sample ID", names(ss), value = TRUE)[1]
if (is.na(colsty) || is.na(colsid)) stop("El sample sheet no tiene columnas de tipo/ID de muestra: ", paste(names(ss), collapse = ", "))
S("Tipos de muestra:"); PT(table(ss[[colsty]], useNA = "ifany"))
# casos: el sample sheet puede traer varios identificadores separados por coma; nos quedamos con el submitter_id TARGET-xx-xxxxxx
ss$case <- sub(",.*$", "", ss[[colcase]])
ss <- ss[grepl("rna_seq", ss[[colfile]]), ]
prim <- grepl("Primary", ss[[colsty]], ignore.case = TRUE) | grepl("-01[A-Z]?(,|$)", ss[[colsid]])
S("Archivos RNA-seq: ", nrow(ss), " | de tumor primario: ", sum(prim))
ss <- ss[prim, ]; ss <- ss[!duplicated(ss$case), ]
S("Casos unicos con RNA-seq de tumor primario: ", nrow(ss))
allf <- unlist(lapply(COUNT_DIRS, function(d) list.files(d, pattern = "star_gene_counts\\.tsv$", recursive = TRUE, full.names = TRUE)))
allf <- allf[!duplicated(basename(allf))]
S("Archivos de conteos encontrados en disco: ", length(allf))
ss$path <- allf[match(ss[[colfile]], basename(allf))]
S("Archivos del sample sheet localizados: ", sum(!is.na(ss$path)), " de ", nrow(ss)); ss <- ss[!is.na(ss$path), ]
rd <- function(p) {
  x <- read.delim(p, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  x <- x[!grepl("^N_", x$gene_id), ]
  tp <- if ("tpm_unstranded" %in% names(x)) x$tpm_unstranded else stop("sin columna tpm_unstranded en ", basename(p))
  tapply(tp, x$gene_name, sum)       # genes duplicados (regiones PAR) se suman
}
L <- lapply(ss$path, rd); genes <- Reduce(intersect, lapply(L, names))
TPM <- sapply(L, function(v) v[genes]); colnames(TPM) <- ss$case
S("Matriz: ", nrow(TPM), " genes x ", ncol(TPM), " muestras")
G <- log2(TPM + 1)
Gf <- G[rowSums(TPM >= 1) >= 0.2 * ncol(TPM), ]
S("Genes expresados (TPM>=1 en >=20% de muestras): ", nrow(Gf))

# -----------------------------------------------------------------------------
SH("3. Union con la clinica y puntuaciones de programas (definicion del script 05)")
gsf <- file.path(RES, "GSE62564_score", "genes_programas.csv")
if (!file.exists(gsf)) stop("Falta ", gsf, " (ejecuta antes el script 05)")
gs <- read.csv(gsf); sets <- split(gs$gen, gs$programa)
sv <- merge(cl, data.frame(case = colnames(Gf)), by = "case"); rownames(sv) <- sv$case
sv <- sv[!is.na(sv$os_days) & !is.na(sv$age_days) & sv$os_days > 0, ]
Zs <- t(scale(t(Gf[, sv$case])))
for (p in names(sets)) { g <- intersect(sets[[p]], rownames(Zs)); sv[[p]] <- colMeans(Zs[g, , drop = FALSE]); S(sprintf("  %-8s genes en el programa %d / presentes %d", p, length(sets[[p]]), length(g))) }
for (p in names(sets)) sv[[paste0(p, "_sd")]] <- as.numeric(scale(sv[[p]]))
sv$logage <- log1p(sv$age_days)
S("Pacientes en el analisis: ", nrow(sv), " | eventos OS: ", sum(sv$os_event), " | EFS: ", sum(sv$efs_event),
  " | MYCN amplificado: ", sum(sv$mycn_amp == 1, na.rm = TRUE), " | estadio 4: ", sum(sv$stage4 == 1, na.rm = TRUE))
write.csv(sv, OUT("puntuaciones_pacientes_TARGET.csv"), row.names = FALSE)
saveRDS(list(Gf = Gf, sv = sv), OUT("cache_TARGET.rds"))   # matriz a nivel gen + clinica: la usa el script 08
S("Comprobacion de sentido biologico (mismas direcciones que en SEQC): puntuacion media segun MYCN y estadio")
for (p in names(sets)) S(sprintf("  %-8s MYCN no amp %.2f, amp %.2f | estadio 4: no %.2f, si %.2f | Wilcoxon p (MYCN) = %.2g", p,
  mean(sv[[p]][sv$mycn_amp == 0], na.rm = TRUE), mean(sv[[p]][sv$mycn_amp == 1], na.rm = TRUE), mean(sv[[p]][sv$stage4 == 0], na.rm = TRUE), mean(sv[[p]][sv$stage4 == 1], na.rm = TRUE),
  tryCatch(wilcox.test(sv[[p]] ~ sv$mycn_amp)$p.value, error = function(e) NA)))

# -----------------------------------------------------------------------------
SH("4. Cox por programa en TARGET-NBL (HR por 1 DE), mismos modelos que en SEQC")
cx1 <- function(form, d, term) {
  m <- tryCatch(coxph(form, data = d), error = function(e) NULL); if (is.null(m)) return(rep(NA_real_, 7))
  s <- summary(m); c(HR = s$conf.int[term, 1], lo = s$conf.int[term, 3], hi = s$conf.int[term, 4], p = s$coefficients[term, 5], n = s$n, ev = s$nevent, C = s$concordance[1]) }
tab <- function(rhs, d = sv, time = "os_days", ev = "os_event") t(sapply(names(sets), function(p) cx1(as.formula(paste0("Surv(", time, ",", ev, ") ~ ", p, "_sd", rhs)), d, paste0(p, "_sd"))))
mods <- list(list("OS univariante", "", sv, "os_days", "os_event"),
             list("OS + edad(log) + MYCN + estadio4", " + logage + mycn_amp + stage4", sv, "os_days", "os_event"),
             list("OS + edad(log) + MYCN + estadio4 + riesgo COG alto", " + logage + mycn_amp + stage4 + HR_cog", sv, "os_days", "os_event"),
             list("OS solo MYCN no amplificado (+ edad + estadio4)", " + logage + stage4", sv[which(sv$mycn_amp == 0), ], "os_days", "os_event"),
             list("EFS univariante", "", sv, "efs_days", "efs_event"),
             list("EFS + edad(log) + MYCN + estadio4", " + logage + mycn_amp + stage4", sv, "efs_days", "efs_event"))
TAB_HR <- list()
for (m in mods) { tt_ <- tab(m[[2]], m[[3]], m[[4]], m[[5]]); S("-- ", m[[1]]); STAB(tt_, 4)
  TAB_HR[[m[[1]]]] <- data.frame(cohorte = "TARGET-NBL", modelo = m[[1]], programa = rownames(tt_), as.data.frame(tt_), row.names = NULL) }
write.csv(do.call(rbind, TAB_HR), OUT("tabla_HR_programas_TARGET.csv"), row.names = FALSE)
S("Referencia SEQC (script 05), CICLO ajustado edad+MYCN+estadio4: OS HR 1.99 (1.57-2.52); solo MYCN no amp 2.04; EFS 1.62.")

SH("5. ?Aporta CICLO sobre el modelo clinico en TARGET-NBL?")
for (ep in list(c("OS", "os_days", "os_event"), c("EFS", "efs_days", "efs_event"))) {
  d <- sv[complete.cases(sv[, c("logage", "mycn_amp", "stage4", "CICLO_sd")]), ]
  m0 <- coxph(as.formula(paste0("Surv(", ep[2], ",", ep[3], ") ~ logage + mycn_amp + stage4")), data = d); m1 <- update(m0, . ~ . + CICLO_sd)
  S(sprintf("%s (n=%d): C clinico = %.3f | C + CICLO = %.3f | LR p = %.3g | HR CICLO por DE = %.2f", ep[1], nrow(d), summary(m0)$concordance[1], summary(m1)$concordance[1], anova(m0, m1)[2, "Pr(>|Chi|)"], exp(coef(m1)["CICLO_sd"])))
}
SH("6. Riesgos proporcionales (cox.zph), OS + edad + MYCN + estadio4")
for (p in names(sets)) {
  z <- tryCatch(cox.zph(coxph(as.formula(paste0("Surv(os_days, os_event) ~ ", p, "_sd + logage + mycn_amp + stage4")), data = sv))$table, error = function(e) NULL)
  if (!is.null(z)) S(sprintf("  %-8s p(programa)=%.3g p(GLOBAL)=%.3g", p, z[1, "p"], z["GLOBAL", "p"]))
}
SH("7. Kaplan-Meier por terciles (OS)")
png(OUT("KM_programas_OS_TARGET.png"), 2200, 1000, res = 200); par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
TAB_S5 <- list()
for (p in names(sets)) {
  ter <- cut(sv[[p]], quantile(sv[[p]], c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("bajo", "medio", "alto"))
  fit <- survfit(Surv(os_days / 365.25, os_event) ~ ter, data = sv); lr <- survdiff(Surv(os_days, os_event) ~ ter, data = sv)
  plot(fit, col = c("#2c7bb6", "grey40", "#d7191c"), lwd = 2, xlab = "anios", ylab = "supervivencia global", main = sprintf("%s (log-rank p=%.2g)", p, 1 - pchisq(lr$chisq, 2)))
  legend("bottomleft", levels(ter), col = c("#2c7bb6", "grey40", "#d7191c"), lwd = 2, bty = "n")
  s5 <- summary(fit, times = 5, extend = TRUE)$surv; S(sprintf("  %-8s S(5 anios) bajo %.2f | medio %.2f | alto %.2f | log-rank p = %.2g", p, s5[1], s5[2], s5[3], 1 - pchisq(lr$chisq, 2)))
  TAB_S5[[p]] <- data.frame(cohorte = "TARGET-NBL", programa = p, tercil = c("bajo", "medio", "alto"), S5 = as.numeric(s5), n = as.integer(table(ter)), logrank_p = 1 - pchisq(lr$chisq, 2))
}
dev.off()
write.csv(do.call(rbind, TAB_S5), OUT("tabla_S5_terciles_TARGET.csv"), row.names = FALSE)

SH("8. Comparacion directa con SEQC dentro de estadio 4 (TARGET-NBL esta enriquecida en estadio 4) y composicion de las cohortes")
st4 <- sv[which(sv$stage4 == 1), ]
S("Pacientes en estadio 4 (TARGET): ", nrow(st4), " de ", nrow(sv), " (", round(100 * nrow(st4) / nrow(sv)), "%) | en SEQC el estadio 4 es el 37%")
S("Muertes (OS): TARGET ", round(100 * mean(sv$os_event)), "% | SEQC 21%")
TAB_ST <- list()
for (p in names(sets)) {
  d <- st4; d$prog_sd <- as.numeric(scale(d[[p]]))
  fnum <- function(rhs, t = "os_days", e = "os_event") { m <- tryCatch(coxph(as.formula(paste0("Surv(", t, ",", e, ") ~ prog_sd", rhs)), data = d), error = function(e) NULL)
    if (is.null(m)) return(c(n = nrow(d), ev = sum(d[[e]]), HR = NA, lo = NA, hi = NA, p = NA)); s <- summary(m)
    c(n = s$n, ev = s$nevent, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) }
  ffmt <- function(v) if (is.na(v["HR"])) "no estimable" else sprintf("n=%d ev=%d HR=%.2f (%.2f-%.2f) p=%.3g", v["n"], v["ev"], v["HR"], v["lo"], v["hi"], v["p"])
  v1 <- fnum(""); v2 <- fnum(" + logage"); v3 <- fnum("", "efs_days", "efs_event")
  S(sprintf("  %-8s estadio 4: OS %s", p, ffmt(v1))); S(sprintf("  %-8s estadio 4: OS + edad(log) %s", "", ffmt(v2))); S(sprintf("  %-8s estadio 4: EFS %s", "", ffmt(v3)))
  TAB_ST[[p]] <- data.frame(cohorte = "TARGET-NBL", programa = p, estrato = "Estadio 4", analisis = c("OS", "OS + edad", "EFS"), rbind(v1, v2, v3), row.names = NULL)
}
write.csv(do.call(rbind, TAB_ST), OUT("tabla_HR_estratos_TARGET.csv"), row.names = FALSE)
S("Referencia SEQC dentro de estadio 4 (script 05, apartado 8): CICLO OS HR 1.50 (1.20-1.89); +edad 1.56 (1.24-1.96); EFS 1.40 (1.15-1.70).")

end_log()
