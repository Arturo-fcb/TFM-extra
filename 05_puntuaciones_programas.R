# =============================================================================
# 05_puntuaciones_programas.R
#   Puntuaciones de PROGRAMAS biologicos (no genes sueltos) en la cohorte SEQC (GSE62564 + GSE49711)
#   Los conjuntos de genes se fijan por ONTOLOGIA (GO, con sus descendientes), no por los resultados de
#   esta cohorte: la pertenencia de cada gen al programa no depende de los datos. Los TERMINOS GO si se
#   eligieron mirando el enriquecimiento (scripts 03/04); por eso la validacion externa es la prueba clave.
#   Programas:  CICLO   = GO:0007059 chromosome segregation
#               RIBO    = GO:0042254 ribosome biogenesis
#               SINAPSIS= GO:0050804 modulation of chemical synaptic transmission
#               INMUNE  = GO:0042110 T cell activation
#   Puntuacion = media del z-score (entre pacientes) del log2 RPM de los genes del programa expresados.
# Salidas  : resultados/GSE62564_score/  (incluye genes_programas.csv, para aplicar la MISMA definicion a otra cohorte)
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("survival", "AnnotationDbi", "org.Hs.eg.db"))
suppressPackageStartupMessages({ library(survival); library(AnnotationDbi); library(org.Hs.eg.db) })
start_log("GSE62564_score")
D <- file.path(DATA, "62564")
num <- function(x) suppressWarnings(as.numeric(x))
PT <- function(x) for (l in capture.output(print(x))) S(l)

# -----------------------------------------------------------------------------
SH("1. Datos: expresion a nivel gen (cache del script 04) y clinica (GSE62564 + GSE49711)")
cache <- file.path(RES, "cache_G_GSE62564.rds")
if (!file.exists(cache)) stop("Falta ", cache, ". Ejecuta antes el script 04 (crea la cache de la matriz a nivel gen).")
G <- readRDS(cache)
meta <- parse_series_matrix(file.path(D, "GSE62564_series_matrix.txt")); meta$id <- sub(" \\[.*\\]$", "", meta$title)
cand <- c(file.path(DATA, "49711", "GSE49711_series_matrix.txt"), file.path(DATA, "49711", "GSE49711_series_matrix.txt.gz"),
          file.path(DATA, "GSE49711_series_matrix.txt"), file.path(DATA, "GSE49711_series_matrix.txt.gz"))
f49 <- cand[file.exists(cand)][1]; if (is.na(f49)) stop("No encuentro GSE49711_series_matrix.txt(.gz)")
m49 <- parse_series_matrix(f49); m49$id <- sub(" \\[.*\\]$", "", m49$title); m49 <- m49[match(meta$id, m49$id), ]
sv <- data.frame(id = meta$id, os_days = num(meta$os.day), os_event = num(meta$os.bin), efs_days = num(meta$efs.day), efs_event = num(meta$efs.bin),
                 age_days = num(meta$age), HR = as.integer(!is.na(meta$high.risk) & meta$high.risk == "HR"), set = meta$training.validation)
sv$mycn_amp <- ifelse(m49$mycn.status %in% c("0", "1"), as.integer(m49$mycn.status), NA)
sv$stage <- factor(trimws(m49$inss.stage), levels = c("1", "2", "3", "4S", "4")); sv$stage4 <- as.integer(sv$stage == "4")
sv <- sv[!is.na(sv$mycn_amp) & !is.na(sv$stage), ]; rownames(sv) <- sv$id; sv$logage <- log1p(sv$age_days)
Gf <- G[rowSums(G >= 0) >= 0.2 * ncol(G), sv$id]
S("Pacientes: ", nrow(sv), " | genes expresados: ", nrow(Gf), " | eventos OS: ", sum(sv$os_event), " | EFS: ", sum(sv$efs_event))

# -----------------------------------------------------------------------------
SH("2. Definicion de los programas (GO, con descendientes) y puntuaciones")
GO <- c(CICLO = "GO:0007059", RIBO = "GO:0042254", SINAPSIS = "GO:0050804", INMUNE = "GO:0042110")
ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unname(GO), keytype = "GOALL", columns = "SYMBOL"))
sets <- lapply(GO, function(go) sort(unique(ann$SYMBOL[ann$GOALL == go & !is.na(ann$SYMBOL)])))
gs_tab <- do.call(rbind, lapply(names(sets), function(n) data.frame(programa = n, GO = GO[[n]], gen = sets[[n]])))
write.csv(gs_tab, OUT("genes_programas.csv"), row.names = FALSE)
Z <- t(scale(t(Gf)))
score <- sapply(sets, function(g) { g <- intersect(g, rownames(Z)); colMeans(Z[g, , drop = FALSE]) })
S("Genes por programa (en GO / expresados en la cohorte): ",
  paste(sprintf("%s %d/%d", names(sets), lengths(sets), sapply(sets, function(g) length(intersect(g, rownames(Z))))), collapse = " | "))
sv <- cbind(sv, score[sv$id, ])
for (p in names(GO)) sv[[paste0(p, "_sd")]] <- as.numeric(scale(sv[[p]]))
write.csv(sv, OUT("puntuaciones_pacientes.csv"), row.names = FALSE)
S("Correlacion Spearman entre programas y con edad:"); PT(round(cor(sv[, c(names(GO), "age_days")], method = "spearman"), 2))

# -----------------------------------------------------------------------------
SH("3. Puntuaciones segun MYCN, estadio y riesgo")
agg <- function(v) tapply(sv[[v]], list(MYCN = sv$mycn_amp), mean)
for (p in names(GO)) {
  S(sprintf("%-8s media | MYCN no amp %.2f, amp %.2f | estadio 4: no %.2f, si %.2f | no-HR %.2f, HR %.2f", p,
            mean(sv[[p]][sv$mycn_amp == 0]), mean(sv[[p]][sv$mycn_amp == 1]), mean(sv[[p]][sv$stage4 == 0]), mean(sv[[p]][sv$stage4 == 1]),
            mean(sv[[p]][sv$HR == 0]), mean(sv[[p]][sv$HR == 1])))
}
S("Wilcoxon p (MYCN amp vs no amp | HR vs no-HR | dentro de MYCN no amp: estadio 4 vs resto):")
for (p in names(GO)) S(sprintf("  %-8s %.2g | %.2g | %.2g", p, wilcox.test(sv[[p]] ~ sv$mycn_amp)$p.value, wilcox.test(sv[[p]] ~ sv$HR)$p.value,
                             wilcox.test(sv[[p]][sv$mycn_amp == 0] ~ sv$stage4[sv$mycn_amp == 0])$p.value))
png(OUT("score_por_grupos.png"), 2200, 1000, res = 200); par(mfrow = c(1, 4), mar = c(7, 4, 3, 1))
grp <- factor(ifelse(sv$mycn_amp == 1, "MYCN amp", ifelse(sv$stage4 == 1, "no amp, E4", "no amp, otros")), levels = c("no amp, otros", "no amp, E4", "MYCN amp"))
for (p in names(GO)) boxplot(sv[[p]] ~ grp, las = 2, main = p, xlab = "", ylab = "puntuacion (z medio)", col = c("grey85", "#9ecae1", "#fc9272"))
dev.off()

# -----------------------------------------------------------------------------
SH("4. Cox por programa (HR por 1 DE de la puntuacion). OS y EFS")
TAB_HR <- list()
cx1 <- function(form, d, term) {
  m <- tryCatch(coxph(form, data = d), error = function(e) NULL); if (is.null(m)) return(rep(NA_real_, 7))
  s <- summary(m); c(HR = s$conf.int[term, 1], lo = s$conf.int[term, 3], hi = s$conf.int[term, 4], p = s$coefficients[term, 5], n = s$n, ev = s$nevent, C = s$concordance[1]) }
tab <- function(rhs, d = sv, time = "os_days", ev = "os_event") t(sapply(names(GO), function(p) cx1(as.formula(paste0("Surv(", time, ",", ev, ") ~ ", p, "_sd", rhs)), d, paste0(p, "_sd"))))
for (m in list(list("OS univariante", "", sv, "os_days", "os_event"),
               list("OS + edad(log) + MYCN + estadio4", " + logage + mycn_amp + stage4", sv, "os_days", "os_event"),
               list("OS + edad(log) + MYCN + estadio4 + HR", " + logage + mycn_amp + stage4 + HR", sv, "os_days", "os_event"),
               list("OS solo MYCN no amplificado (+ edad + estadio4)", " + logage + stage4", sv[sv$mycn_amp == 0, ], "os_days", "os_event"),
               list("OS solo TRAINING (+ edad + MYCN + estadio4)", " + logage + mycn_amp + stage4", sv[sv$set == "Training", ], "os_days", "os_event"),
               list("OS solo VALIDATION (+ edad + MYCN + estadio4)", " + logage + mycn_amp + stage4", sv[sv$set == "Validation", ], "os_days", "os_event"),
               list("EFS univariante", "", sv, "efs_days", "efs_event"),
               list("EFS + edad(log) + MYCN + estadio4", " + logage + mycn_amp + stage4", sv, "efs_days", "efs_event"))) {
  tt_ <- tab(m[[2]], m[[3]], m[[4]], m[[5]]); S("-- ", m[[1]]); STAB(tt_, 4)
  TAB_HR[[m[[1]]]] <- data.frame(cohorte = "SEQC", modelo = m[[1]], programa = rownames(tt_), as.data.frame(tt_), row.names = NULL)
}
write.csv(do.call(rbind, TAB_HR), OUT("tabla_HR_programas_SEQC.csv"), row.names = FALSE)
S("Todos los programas a la vez (OS, + edad + MYCN + estadio4): ")
m_all <- coxph(Surv(os_days, os_event) ~ CICLO_sd + RIBO_sd + SINAPSIS_sd + INMUNE_sd + logage + mycn_amp + stage4, data = sv)
STAB(round(summary(m_all)$coefficients[, c(2, 5)], 4), 4)

# -----------------------------------------------------------------------------
SH("5. ?Aporta el programa de ciclo celular informacion sobre el modelo clinico (edad + MYCN + estadio)?")
for (ep in list(c("OS", "os_days", "os_event"), c("EFS", "efs_days", "efs_event"))) {
  f0 <- as.formula(paste0("Surv(", ep[2], ",", ep[3], ") ~ logage + mycn_amp + stage4"))
  m0 <- coxph(f0, data = sv); m1 <- update(m0, . ~ . + CICLO_sd)
  S(sprintf("%s: C clinico = %.3f | C + CICLO = %.3f | LR p = %.3g | coef CICLO (por DE) HR = %.2f", ep[1], summary(m0)$concordance[1], summary(m1)$concordance[1],
            anova(m0, m1)[2, "Pr(>|Chi|)"], exp(coef(m1)["CICLO_sd"])))
}

# -----------------------------------------------------------------------------
SH("6. Riesgos proporcionales (cox.zph) del modelo OS + edad + MYCN + estadio4, por programa")
for (p in names(GO)) {
  z <- cox.zph(coxph(as.formula(paste0("Surv(os_days, os_event) ~ ", p, "_sd + logage + mycn_amp + stage4")), data = sv))$table
  S(sprintf("  %-8s p(programa)=%.3g p(edad)=%.3g p(MYCN)=%.3g p(estadio4)=%.3g p(GLOBAL)=%.3g", p, z[1, "p"], z["logage", "p"], z["mycn_amp", "p"], z["stage4", "p"], z["GLOBAL", "p"]))
}

# -----------------------------------------------------------------------------
SH("7. Kaplan-Meier por terciles de cada programa (OS): log-rank y supervivencia a 5 anios")
png(OUT("KM_programas_OS.png"), 2200, 1000, res = 200); par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
TAB_S5 <- list()
for (p in names(GO)) {
  ter <- cut(sv[[p]], quantile(sv[[p]], c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("bajo", "medio", "alto"))
  fit <- survfit(Surv(os_days / 365.25, os_event) ~ ter, data = sv); lr <- survdiff(Surv(os_days, os_event) ~ ter, data = sv)
  plot(fit, col = c("#2c7bb6", "grey40", "#d7191c"), lwd = 2, xlab = "anios", ylab = "supervivencia global", main = sprintf("%s (log-rank p=%.2g)", p, 1 - pchisq(lr$chisq, 2)))
  legend("bottomleft", levels(ter), col = c("#2c7bb6", "grey40", "#d7191c"), lwd = 2, bty = "n")
  s5 <- summary(fit, times = 5, extend = TRUE)$surv; S(sprintf("  %-8s S(5 anios) bajo %.2f | medio %.2f | alto %.2f | log-rank p = %.2g", p, s5[1], s5[2], s5[3], 1 - pchisq(lr$chisq, 2)))
  TAB_S5[[p]] <- data.frame(cohorte = "SEQC", programa = p, tercil = c("bajo", "medio", "alto"), S5 = as.numeric(s5), n = as.integer(table(ter)), logrank_p = 1 - pchisq(lr$chisq, 2))
}
dev.off()
write.csv(do.call(rbind, TAB_S5), OUT("tabla_S5_terciles_SEQC.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
SH("8. Efecto del programa DENTRO de estratos clinicos (la puntuacion se re-estandariza dentro de cada estrato)")
st_num <- function(d, rhs, ep = c("os_days", "os_event")) {
  d$prog_sd <- as.numeric(scale(d$prog))
  m <- tryCatch(coxph(as.formula(paste0("Surv(", ep[1], ",", ep[2], ") ~ prog_sd", rhs)), data = d), error = function(e) NULL)
  if (is.null(m)) return(c(n = nrow(d), ev = sum(d[[ep[2]]]), HR = NA, lo = NA, hi = NA, p = NA)); s <- summary(m)
  c(n = s$n, ev = s$nevent, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) }
st_fmt <- function(v) if (is.na(v["HR"])) "no estimable" else sprintf("n=%d ev=%d HR=%.2f (%.2f-%.2f) p=%.3g", v["n"], v["ev"], v["HR"], v["lo"], v["hi"], v["p"])
st_cox <- function(d, rhs, ep = c("os_days", "os_event")) st_fmt(st_num(d, rhs, ep))
TAB_ST <- list()
strata <- list("Todos" = sv, "Estadio 4" = sv[sv$stage4 == 1, ], "Alto riesgo (HR)" = sv[sv$HR == 1, ], "No alto riesgo" = sv[sv$HR == 0, ],
               "Estadio 4 y MYCN no amplificado" = sv[sv$stage4 == 1 & sv$mycn_amp == 0, ], "HR y MYCN amplificado" = sv[sv$HR == 1 & sv$mycn_amp == 1, ])
for (p in c("CICLO", "RIBO")) {
  S("-- Programa ", p, " (HR por 1 DE)")
  for (n in names(strata)) { d <- strata[[n]]; d$prog <- d[[p]]
    v1 <- st_num(d, ""); v2 <- st_num(d, " + logage"); v3 <- st_num(d, "", c("efs_days", "efs_event"))
    S(sprintf("  %-34s OS: %s", n, st_fmt(v1))); S(sprintf("  %-34s OS + edad(log): %s", "", st_fmt(v2)))
    S(sprintf("  %-34s EFS: %s", "", st_fmt(v3)))
    TAB_ST[[paste(p, n)]] <- data.frame(cohorte = "SEQC", programa = p, estrato = n, analisis = c("OS", "OS + edad", "EFS"), rbind(v1, v2, v3), row.names = NULL) }
}
write.csv(do.call(rbind, TAB_ST), OUT("tabla_HR_estratos_SEQC.csv"), row.names = FALSE)
S("CICLO medio por estatus MYCN dentro de estadio 4: no amp ", round(mean(sv$CICLO[sv$stage4 == 1 & sv$mycn_amp == 0]), 2), " | amp ", round(mean(sv$CICLO[sv$stage4 == 1 & sv$mycn_amp == 1]), 2))
S("Composicion de SEQC: estadio 4 = ", round(100 * mean(sv$stage4 == 1)), "% | alto riesgo = ", round(100 * mean(sv$HR == 1)), "% | muertes = ", round(100 * mean(sv$os_event)), "%")

end_log()
