# =============================================================================
# 07_figuras_finales.R   FIGURAS FINALES EN CASTELLANO (para la memoria y el correo al tutor)
#   Todas las figuras salen de los CSV que dejan los scripts 01-06 (nada se escribe a mano):
#     01-04 -> DE_*.csv, Cox_*.csv, RMA_expr.rds, *_ORA_GOBP_*_simplificado.csv
#     05    -> tabla_HR_programas_SEQC.csv, tabla_HR_estratos_SEQC.csv, tabla_S5_terciles_SEQC.csv
#     06    -> tabla_HR_programas_TARGET.csv, tabla_HR_estratos_TARGET.csv, tabla_S5_terciles_TARGET.csv
#   Si falta una entrada, esa figura se omite con un aviso (el resto sigue).
# Salidas  : resultados/figuras/*.png  (+ tablas .csv de apoyo y RESUMEN_figuras.txt)
# Requiere : ggplot2, ggrepel, scales
# Los terminos GO se traducen con scripts/traduccion_GO.csv (columnas: term, es). Los que falten salen en
# ingles y quedan listados en el RESUMEN para completarlos.
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("ggplot2", "ggrepel", "scales"))
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
start_log("figuras")
PT <- function(x) for (l in capture.output(print(x))) S(l)

# ---- estilo comun ----------------------------------------------------------
COL_TARGET <- "#2C7BB6"; COL_UP <- "#D7191C"; COL_DOWN <- "#2C7BB6"
COLS <- c("SEQC" = "grey55", "TARGET-NBL" = COL_TARGET)
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")
tema <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10.5), legend.position = "bottom",
        plot.title.position = "plot", plot.caption = element_text(colour = "grey40", size = 9, hjust = 0))
guardar <- function(p, nombre, w = 9, h = 5.5) {
  f <- OUT(nombre)
  if (has("ragg")) ggsave(f, p, width = w, height = h, dpi = 200, bg = "white", device = ragg::agg_png) else ggsave(f, p, width = w, height = h, dpi = 200, bg = "white")
  S("  figura: ", nombre) }
leer <- function(...) { f <- file.path(RES, ...); if (!file.exists(f)) { S("  (falta ", f, ")"); return(NULL) }; read.csv(f, check.names = FALSE) }
omitida <- function(nombre, motivo) S("  OMITIDA ", nombre, ": ", motivo)

# -----------------------------------------------------------------------------
SH("1. Concordancia entre contrastes de los tres datasets (genes comunes, Spearman de t)")
gn <- leer("GSE18271", "DE_GNB_vs_GN_gen.csv"); cl <- leer("GSE212877", "D_pareado_mixto.csv")
hr <- leer("GSE62564", "DE_HR_vs_nonHR_gen.csv"); cx <- leer("GSE62564", "Cox_transcriptoma_completo.csv")
if (!any(sapply(list(gn, cl, hr, cx), is.null))) {
  tv <- function(d) setNames(d[["t"]], d[[1]])
  gn_t <- tv(gn); cl_t <- tv(cl); hr_t <- tv(hr)
  cx_z <- setNames(sign(cx$coef) * (-qnorm(cx$p / 2)), cx[[1]])
  comunes <- Reduce(intersect, list(names(gn_t), names(cl_t), names(hr_t)))
  rho <- function(a, b) cor(a[comunes], b[comunes], method = "spearman")
  conc <- data.frame(
    comparacion = c("GN/GNB frente a l\u00edneas celulares", "GN/GNB frente a alto riesgo", "L\u00edneas celulares frente a alto riesgo", "Alto riesgo frente a Cox (mismo dataset)"),
    grupo = c(rep("Entre datasets distintos", 3), "Dentro de SEQC (control)"),
    rho = c(rho(gn_t, cl_t), rho(gn_t, hr_t), rho(cl_t, hr_t), rho(hr_t, cx_z)))
  S("Genes comunes a los tres datasets: ", length(comunes)); PT(round(conc$rho, 3))
  write.csv(conc, OUT("tabla_concordancia_datasets.csv"), row.names = FALSE)
  conc$comparacion <- factor(conc$comparacion, levels = rev(conc$comparacion))
  p <- ggplot(conc, aes(rho, comparacion, fill = grupo)) + geom_col(width = .55) +
    geom_vline(xintercept = 0, colour = "grey30") +
    geom_text(aes(label = fmt(rho), hjust = ifelse(rho < 0, 1.15, -0.15)), size = 3.8) +
    scale_fill_manual(values = c("Entre datasets distintos" = COL_TARGET, "Dentro de SEQC (control)" = "grey65"), name = NULL) +
    scale_x_continuous(limits = c(-0.6, 1), labels = function(x) sub("-", "\u2212", fmt(x, 1))) +
    labs(title = "Los datasets casi no se parecen entre s\u00ed; solo coinciden dos medidas de SEQC",
         subtitle = sprintf("Correlaci\u00f3n de Spearman de los estad\u00edsticos t entre contrastes (%s genes comunes)", format(length(comunes), big.mark = ".")),
         x = "\u03c1 de Spearman", y = NULL) + tema
  guardar(p, "fig_concordancia_datasets.png", 9, 4.2)
} else omitida("fig_concordancia_datasets.png", "faltan CSV de los scripts 01-03")

# -----------------------------------------------------------------------------
SH("2. Efecto del programa CICLO: por estratos clinicos y SEQC frente a TARGET-NBL")
hs <- leer("GSE62564_score", "tabla_HR_estratos_SEQC.csv"); ht <- leer("TARGET_NBL", "tabla_HR_estratos_TARGET.csv")
ps <- leer("GSE62564_score", "tabla_HR_programas_SEQC.csv"); pt <- leer("TARGET_NBL", "tabla_HR_programas_TARGET.csv")
forest <- function(d, titulo, subtitulo, nombre, h = 5, ticks = c(1, 1.5, 2, 3)) {
  d$fila <- factor(d$fila, levels = rev(unique(d$fila)))
  d$etq <- sprintf("%s (%s\u2013%s)", fmt(d$HR), fmt(d$lo), fmt(d$hi))
  pos <- if (length(unique(d$cohorte)) > 1) position_dodge(width = .6) else position_identity()
  p <- ggplot(d, aes(HR, fila, colour = cohorte, group = cohorte)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.9, position = pos) +
    geom_point(size = 3, position = pos) +
    geom_text(aes(x = hi, label = etq), hjust = -0.12, size = 3.2, colour = "black", position = pos) +
    scale_colour_manual(values = COLS, name = NULL) +
    scale_x_log10(breaks = ticks, labels = function(x) sub("\\.", ",", format(x)), expand = expansion(mult = c(0.03, 0.45))) +
    labs(title = titulo, subtitle = subtitulo, x = "HR por desviaci\u00f3n est\u00e1ndar de CICLO (escala logar\u00edtmica); l\u00ednea discontinua: HR = 1", y = NULL) + tema
  guardar(p, nombre, 9.5, h)
}
sel <- function(d, estr, an = "OS", prog = "CICLO") d[d$programa == prog & d$estrato == estr & d$analisis == an, ]
if (!is.null(hs)) {
  filas <- list(list("SEQC, todos los pacientes", sel(hs, "Todos")), list("SEQC, alto riesgo", sel(hs, "Alto riesgo (HR)")),
                list("SEQC, estadio 4", sel(hs, "Estadio 4")), list("SEQC, estadio 4 sin MYCN amplificado", sel(hs, "Estadio 4 y MYCN no amplificado")))
  if (!is.null(ht)) filas[[5]] <- list("TARGET-NBL, estadio 4", sel(ht, "Estadio 4"))
  d2 <- do.call(rbind, lapply(filas, function(f) cbind(fila = sprintf("%s\n(n = %d, %d muertes)", f[[1]], f[[2]]$n, f[[2]]$ev), f[[2]][, c("cohorte", "HR", "lo", "hi")])))
  PT(d2); write.csv(d2, OUT("tabla_forest_CICLO_estratos.csv"), row.names = FALSE)
  forest(d2, "El efecto de CICLO se debilita al acotar el riesgo cl\u00ednico y desaparece en TARGET-NBL",
         "HR por desviaci\u00f3n est\u00e1ndar de la puntuaci\u00f3n CICLO, supervivencia global sin ajustar, IC 95 %", "fig_forest_CICLO_estratos.png", 5.6)
} else omitida("fig_forest_CICLO_estratos.png", "falta tabla_HR_estratos_SEQC.csv (ejecuta el script 05)")

if (!is.null(hs) && !is.null(ht) && !is.null(ps) && !is.null(pt)) {
  mod <- "OS + edad(log) + MYCN + estadio4"
  g <- function(d, modelo) d[d$programa == "CICLO" & d$modelo == modelo, c("cohorte", "HR", "lo", "hi", "n", "ev")]
  filas <- list(
    list("SG, cohorte completa\najustada por edad, MYCN y estadio 4", rbind(g(ps, mod), g(pt, mod))),
    list("SG, solo estadio 4\nsin ajustar", rbind(cbind(sel(hs, "Estadio 4")[, c("cohorte", "HR", "lo", "hi", "n")], ev = sel(hs, "Estadio 4")$ev), cbind(sel(ht, "Estadio 4")[, c("cohorte", "HR", "lo", "hi", "n")], ev = sel(ht, "Estadio 4")$ev))),
    list("SG, solo estadio 4\najustada por edad", rbind(cbind(sel(hs, "Estadio 4", "OS + edad")[, c("cohorte", "HR", "lo", "hi", "n")], ev = sel(hs, "Estadio 4", "OS + edad")$ev), cbind(sel(ht, "Estadio 4", "OS + edad")[, c("cohorte", "HR", "lo", "hi", "n")], ev = sel(ht, "Estadio 4", "OS + edad")$ev))),
    list("SLE, solo estadio 4\nsin ajustar", rbind(cbind(sel(hs, "Estadio 4", "EFS")[, c("cohorte", "HR", "lo", "hi", "n")], ev = sel(hs, "Estadio 4", "EFS")$ev), cbind(sel(ht, "Estadio 4", "EFS")[, c("cohorte", "HR", "lo", "hi", "n")], ev = sel(ht, "Estadio 4", "EFS")$ev))))
  d3 <- do.call(rbind, lapply(filas, function(f) cbind(fila = f[[1]], f[[2]][, c("cohorte", "HR", "lo", "hi")])))
  PT(d3); write.csv(d3, OUT("tabla_forest_CICLO_SEQC_vs_TARGET.csv"), row.names = FALSE)
  forest(d3, "El efecto de CICLO en SEQC desaparece en TARGET-NBL",
         "HR por desviaci\u00f3n est\u00e1ndar de CICLO, IC 95 %. SG: supervivencia global; SLE: supervivencia libre de evento", "fig_forest_CICLO_SEQC_vs_TARGET.png", 5.8)

  # todos los programas, OS ajustada (comparacion SEQC / TARGET)
  dp <- rbind(ps[ps$modelo == mod, ], pt[pt$modelo == mod, ])
  dp$fila <- dp$programa; dp <- dp[, c("fila", "cohorte", "HR", "lo", "hi")]
  forest(dp, "Cuatro programas: ninguno se asocia con la supervivencia en TARGET-NBL tras ajustar",
         "HR por desviaci\u00f3n est\u00e1ndar, supervivencia global ajustada por edad, MYCN y estadio 4, IC 95 %", "fig_forest_4programas_SEQC_vs_TARGET.png", 5, ticks = c(0.5, 1, 2))
} else omitida("fig_forest_CICLO_SEQC_vs_TARGET.png", "faltan las tablas de los scripts 05 y 06")

# -----------------------------------------------------------------------------
SH("3. Supervivencia a 5 anios por terciles de cada programa (Kaplan-Meier, sin ajustar)")
s5s <- leer("GSE62564_score", "tabla_S5_terciles_SEQC.csv"); s5t <- leer("TARGET_NBL", "tabla_S5_terciles_TARGET.csv")
if (!is.null(s5s) && !is.null(s5t)) {
  d <- rbind(s5s, s5t); d$tercil <- factor(d$tercil, levels = c("bajo", "medio", "alto"), labels = c("Tercil bajo", "Tercil medio", "Tercil alto"))
  d$cohorte <- factor(d$cohorte, levels = c("SEQC", "TARGET-NBL")); PT(d)
  write.csv(d, OUT("tabla_S5_terciles_SEQC_TARGET.csv"), row.names = FALSE)
  dd <- d[d$programa == "CICLO", ]
  p <- ggplot(dd, aes(tercil, S5, colour = cohorte, group = cohorte)) + geom_line(linewidth = 1) + geom_point(size = 3.2) +
    geom_text(aes(label = scales::percent(S5, accuracy = 1, suffix = " %"), vjust = ifelse(cohorte == "SEQC", -1.1, 2.0)), size = 3.6, colour = "black") +
    scale_colour_manual(values = COLS, name = NULL) + scale_y_continuous(labels = function(x) paste0(round(100 * x), " %"), limits = c(0, 1.02)) +
    labs(title = "La supervivencia cae con el ciclo celular en SEQC, no en TARGET-NBL",
         subtitle = "Supervivencia global a 5 a\u00f1os por terciles de la puntuaci\u00f3n CICLO dentro de cada cohorte (Kaplan-Meier sin ajustar)",
         x = "de menos a m\u00e1s puntuaci\u00f3n de ciclo celular", y = "Supervivencia global a 5 a\u00f1os") + tema
  guardar(p, "fig_supervivencia_terciles_CICLO.png", 8.5, 5)
  p4 <- ggplot(d, aes(tercil, S5, colour = cohorte, group = cohorte)) + geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
    facet_wrap(~ programa, nrow = 1) + scale_colour_manual(values = COLS, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), " %"), limits = c(0, 1.02)) +
    labs(title = "Supervivencia a 5 a\u00f1os por terciles de cada programa", subtitle = "Kaplan-Meier sin ajustar, dentro de cada cohorte", x = NULL, y = "Supervivencia global a 5 a\u00f1os") +
    tema + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  guardar(p4, "fig_supervivencia_terciles_4programas.png", 11, 4.6)
} else omitida("fig_supervivencia_terciles_CICLO.png", "faltan tabla_S5_terciles_* (scripts 05 y 06)")

# -----------------------------------------------------------------------------
SH("4. PCA de GSE18271 con las 22 muestras etiquetadas")
rds <- file.path(RES, "GSE18271", "RMA_expr.rds"); sm <- file.path(DATA, "18271", "GSE18271_series_matrix.txt")
if (file.exists(rds) && file.exists(sm)) {
  expr <- readRDS(rds); meta <- parse_series_matrix(sm)
  grp <- factor(ifelse(grepl("blastoma", meta$tissue), "Ganglioneuroblastoma (GNB)", "Ganglioneuroma (GN)"), levels = c("Ganglioneuroma (GN)", "Ganglioneuroblastoma (GNB)"))
  stopifnot(ncol(expr) == nrow(meta))
  pc <- prcomp(t(expr), scale. = TRUE); ve <- summary(pc)$importance[2, ] * 100
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], grupo = grp, muestra = meta$geo_accession)
  S("Muestras en el PCA: ", nrow(d), " | PC1 ", round(ve[1], 1), "% | PC2 ", round(ve[2], 1), "%")
  p <- ggplot(d, aes(PC1, PC2, colour = grupo, label = muestra)) + geom_point(size = 3.5) +
    geom_text_repel(size = 3.1, max.overlaps = Inf, show.legend = FALSE, min.segment.length = 0.1) +
    scale_colour_manual(values = c("Ganglioneuroma (GN)" = "#1B9E77", "Ganglioneuroblastoma (GNB)" = "#D95F02"), name = NULL) +
    labs(title = "PCA de GSE18271: GN y GNB se separan en parte, con solapamiento",
         subtitle = sprintf("22 muestras, todas las sondas tras RMA. PC1 %s %%, PC2 %s %%", fmt(ve[1], 1), fmt(ve[2], 1)),
         x = sprintf("PC1 (%s %%)", fmt(ve[1], 1)), y = sprintf("PC2 (%s %%)", fmt(ve[2], 1))) + tema
  guardar(p, "fig_PCA_GSE18271_22_muestras.png", 8.5, 6.2)
} else omitida("fig_PCA_GSE18271_22_muestras.png", "falta RMA_expr.rds (script 01) o la series matrix de GSE18271")

# -----------------------------------------------------------------------------
SH("5. Volcanos con umbrales visibles")
volcan <- function(tt, titulo, nombre, grupo_up, grupo_down, xlab, nota_extremos = FALSE) {
  d <- tt; d$gen <- d[[1]]; d$ml <- -log10(pmax(d$adj.P.Val, 1e-300))
  d$cat <- "Sin cambio relevante"
  d$cat[d$adj.P.Val < FDR_CUT & d$logFC >  LFC_CUT] <- paste("Mayor en", grupo_up)
  d$cat[d$adj.P.Val < FDR_CUT & d$logFC < -LFC_CUT] <- paste("Mayor en", grupo_down)
  n_up <- sum(d$cat == paste("Mayor en", grupo_up)); n_dw <- sum(d$cat == paste("Mayor en", grupo_down))
  cols <- setNames(c(COL_UP, COL_DOWN, "grey75"), c(paste("Mayor en", grupo_up), paste("Mayor en", grupo_down), "Sin cambio relevante"))
  ymax <- max(d$ml); xmax <- max(abs(d$logFC))
  sub <- sprintf("%s genes testados; umbrales: FDR < %s y |log2FC| > %d. %s genes con mayor expresi\u00f3n en %s y %s en %s",
                 format(nrow(d), big.mark = "."), fmt(FDR_CUT, 2), LFC_CUT, format(n_up, big.mark = "."), grupo_up, format(n_dw, big.mark = "."), grupo_down)
  if (nota_extremos) {
    ext <- d[abs(d$logFC) > 8, ]
    sub <- paste0(sub, sprintf("\nColas con |log2FC| > 8: %d genes; el %d %% tiene expresi\u00f3n media < 1 log2 CPM (log2FC inflado por conteos muy bajos en un grupo)",
                               nrow(ext), round(100 * mean(ext$AveExpr < 1))))
    S("  ", nombre, ": ", nrow(ext), " genes con |log2FC|>8; % con AveExpr<1 = ", round(100 * mean(ext$AveExpr < 1)))
  }
  p <- ggplot(d, aes(logFC, ml, colour = cat)) + geom_point(alpha = .5, size = 1) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey30") +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", colour = "grey30") +
    annotate("text", x = xmax * 0.98, y = -log10(FDR_CUT), label = "FDR = 0,05", vjust = -0.6, hjust = 1, size = 3.3, colour = "grey20") +
    annotate("text", x = LFC_CUT, y = ymax * 0.98, label = "log2FC = 1", hjust = -0.1, size = 3.3, colour = "grey20") +
    annotate("text", x = -LFC_CUT, y = ymax * 0.98, label = "log2FC = \u22121", hjust = 1.1, size = 3.3, colour = "grey20") +
    scale_colour_manual(values = cols, name = NULL) + labs(title = titulo, subtitle = sub, x = xlab, y = "\u2212log10(FDR)") + tema
  guardar(p, nombre, 9, 6.2)
}
d1 <- leer("GSE18271", "DE_GNB_vs_GN_gen.csv")
if (!is.null(d1)) volcan(d1, "GNB frente a GN (GSE18271, nivel de gen)", "fig_volcano_GSE18271.png", "GNB", "GN", "log2FC (GNB / GN)")
d2 <- leer("GSE212877", "A_SHEPvsADRN4_voom.csv")
if (!is.null(d2)) volcan(d2, "SH-EP frente a l\u00edneas adren\u00e9rgicas (GSE212877, contraste A)", "fig_volcano_GSE212877_A.png", "SH-EP", "ADRN", "log2FC (SH-EP / ADRN)", TRUE)
d3 <- leer("GSE212877", "D_pareado_mixto.csv")
if (!is.null(d3)) volcan(d3, "L\u00edneas S frente a N, pareado con l\u00ednea como efecto aleatorio (GSE212877)", "fig_volcano_GSE212877_D.png", "l\u00edneas S", "l\u00edneas N", "log2FC (S / N)", TRUE)
d4 <- leer("GSE62564", "DE_HR_vs_nonHR_gen.csv")
if (!is.null(d4)) volcan(d4, "Alto riesgo frente a no alto riesgo (GSE62564)", "fig_volcano_GSE62564_HR.png", "alto riesgo", "no alto riesgo", "log2FC (alto riesgo / no alto riesgo)")
d5 <- leer("GSE62564_ajuste", "DE_HR_ajustado_MYCN.csv")
if (!is.null(d5)) volcan(d5, "Alto riesgo frente a no alto riesgo, ajustado por MYCN (GSE62564)", "fig_volcano_GSE62564_HR_ajustado_MYCN.png", "alto riesgo", "no alto riesgo", "log2FC (alto riesgo / no alto riesgo)")

# -----------------------------------------------------------------------------
SH("6. Graficos de puntos de GO (ORA) en castellano")
trad_f <- file.path(ROOT, "scripts", "traduccion_GO.csv")
trad <- if (file.exists(trad_f)) read.csv(trad_f, stringsAsFactors = FALSE, fileEncoding = "UTF-8") else data.frame(term = character(), es = character())
faltan <- character()
es <- function(x) { r <- trad$es[match(x, trad$term)]; faltan <<- union(faltan, x[is.na(r)]); ifelse(is.na(r), x, r) }
wrap <- function(x, w = 46) vapply(x, function(s) paste(strwrap(s, w), collapse = "\n"), "")
gof <- function(csv_dir, csv, titulo, sub, nombre, n = 12) {
  d <- leer(csv_dir, csv); if (is.null(d) || !nrow(d)) { omitida(nombre, "falta el CSV de entrada"); return(invisible()) }
  d <- head(d[order(d$p.adjust), ], n)
  d$razon <- vapply(strsplit(d$GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), 1)
  d$termino <- wrap(es(d$Description)); d$termino <- factor(d$termino, levels = rev(unique(d$termino)))
  p <- ggplot(d, aes(razon, termino, size = Count, colour = -log10(p.adjust))) + geom_point() +
    scale_colour_viridis_c(option = "D", direction = -1, name = "\u2212log10(FDR)") + scale_size_continuous(name = "N\u00ba de genes", range = c(3, 8)) +
    scale_x_continuous(labels = function(x) fmt(x, 2)) +
    labs(title = titulo, subtitle = sub, x = "Proporci\u00f3n de genes del conjunto presentes en la lista", y = NULL,
         caption = "ORA con Gene Ontology (procesos biol\u00f3gicos, t\u00e9rminos simplificados); FDR por Benjamini-Hochberg") + tema + theme(legend.position = "right")
  guardar(p, nombre, 9.5, 5.8)
}
gof("GSE18271", "GNBvsGN_ORA_GOBP_up_simplificado.csv", "GNB frente a GN: procesos de los genes m\u00e1s expresados en GNB", "GSE18271, ORA de los genes con FDR < 0,05 y log2FC > 1", "fig_GO_GSE18271_mas_en_GNB.png")
gof("GSE18271", "GNBvsGN_ORA_GOBP_down_simplificado.csv", "GNB frente a GN: procesos de los genes m\u00e1s expresados en GN", "GSE18271, ORA de los genes con FDR < 0,05 y log2FC < \u22121", "fig_GO_GSE18271_mas_en_GN.png")
gof("GSE212877", "D_pareado_mixto_ORA_GOBP_up_simplificado.csv", "L\u00edneas S frente a N: procesos de los genes m\u00e1s expresados en S", "GSE212877, contraste pareado con l\u00ednea como efecto aleatorio", "fig_GO_GSE212877_mas_en_S.png")
gof("GSE212877", "D_pareado_mixto_ORA_GOBP_down_simplificado.csv", "L\u00edneas S frente a N: procesos de los genes m\u00e1s expresados en N", "GSE212877, contraste pareado con l\u00ednea como efecto aleatorio", "fig_GO_GSE212877_mas_en_N.png")
gof("GSE62564_ajuste", "HR_ajustado_MYCN_ORA_GOBP_up_simplificado.csv", "Alto riesgo ajustado por MYCN: procesos de los genes m\u00e1s expresados", "GSE62564, alto riesgo frente a no alto riesgo ajustando por estatus de MYCN", "fig_GO_GSE62564_HR_ajustado_mas_en_HR.png")
gof("GSE62564_ajuste", "HR_ajustado_MYCN_ORA_GOBP_down_simplificado.csv", "Alto riesgo ajustado por MYCN: procesos de los genes menos expresados", "GSE62564, alto riesgo frente a no alto riesgo ajustando por estatus de MYCN", "fig_GO_GSE62564_HR_ajustado_menos_en_HR.png")
gof("GSE62564_ajuste", "M1_edad_MYCN_est_ORA_GOBP_up_simplificado.csv", "Genes asociados a peor supervivencia (Cox ajustado): procesos", "GSE62564, genes con FDR < 0,05 y coeficiente positivo, ajustado por edad, MYCN y estadio 4", "fig_GO_GSE62564_Cox_riesgo.png")
if (length(faltan)) { S("TERMINOS GO SIN TRADUCCION (salen en ingles; anadelos a scripts/traduccion_GO.csv):"); for (t in faltan) S("   ", t) } else S("Todos los terminos GO mostrados estan traducidos.")

# -----------------------------------------------------------------------------
SH("7. Kaplan-Meier de los cuatro genes candidatos (corte en la mediana, supervivencia global, n = 498)")
cache <- file.path(RES, "cache_G_GSE62564.rds"); f62 <- file.path(DATA, "62564", "GSE62564_series_matrix.txt")
if (file.exists(cache) && file.exists(f62) && requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  G <- readRDS(cache); num <- function(x) suppressWarnings(as.numeric(x))
  mt <- parse_series_matrix(f62); mt$id <- sub(" \\[.*\\]$", "", mt$title)
  sv <- data.frame(id = mt$id, t = num(mt$os.day) / 365.25, e = num(mt$os.bin)); sv <- sv[sv$id %in% colnames(G), ]
  genes <- c("ARPP21", "DLG2", "AXL", "MYCN"); curvas <- list(); tab <- list()
  for (g in genes) {
    d <- sv; d$x <- as.numeric(G[g, d$id]); d$grp <- factor(ifelse(d$x >= median(d$x), "Alto", "Bajo"), levels = c("Bajo", "Alto"))
    fit <- survfit(Surv(t, e) ~ grp, data = d); lr <- survdiff(Surv(t, e) ~ grp, data = d); p <- 1 - pchisq(lr$chisq, 1)
    sm <- summary(fit, times = c(5, 10), extend = TRUE)
    tab[[g]] <- data.frame(gen = g, n_bajo = sum(d$grp == "Bajo"), n_alto = sum(d$grp == "Alto"),
                           S10_bajo = sm$surv[sm$strata == "grp=Bajo" & sm$time == 10], S10_alto = sm$surv[sm$strata == "grp=Alto" & sm$time == 10], logrank_p = p)
    pl <- sprintf("%s (log-rank p = %s)", g, if (p < 1e-4) formatC(p, format = "e", digits = 1, decimal.mark = ",") else formatC(p, format = "f", digits = 4, decimal.mark = ","))
    st <- rep(names(fit$strata), fit$strata)
    curvas[[g]] <- rbind(data.frame(panel = pl, grupo = "Bajo", tiempo = 0, s = 1)[0, ],
      data.frame(panel = pl, grupo = sub("grp=", "", st), tiempo = fit$time, s = fit$surv, lo = fit$lower, hi = fit$upper)[, c("panel", "grupo", "tiempo", "s")])
    curvas[[g]] <- rbind(data.frame(panel = pl, grupo = c("Bajo", "Alto"), tiempo = 0, s = 1), curvas[[g]])
  }
  tab <- do.call(rbind, tab); PT(round(tab[, -1], 4)); write.csv(tab, OUT("tabla_KM_4genes.csv"), row.names = FALSE)
  dd <- do.call(rbind, curvas); dd$panel <- factor(dd$panel, levels = unique(dd$panel)); dd$grupo <- factor(dd$grupo, levels = c("Alto", "Bajo"), labels = c("Expresión alta", "Expresión baja"))
  p <- ggplot(dd, aes(tiempo, s, colour = grupo)) + geom_step(linewidth = 0.9) + facet_wrap(~ panel, nrow = 2) +
    scale_colour_manual(values = c("Expresión alta" = "#D7191C", "Expresión baja" = "#2C7BB6"), name = NULL) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), " %")) +
    labs(title = "Supervivencia global según la expresión de cada gen candidato", subtitle = "GSE62564 (n = 498), corte en la mediana del gen; Kaplan-Meier sin ajustar",
         x = "Tiempo (años)", y = "Supervivencia global") + tema + theme(strip.text = element_text(face = "bold.italic", size = 11))
  guardar(p, "fig_KM_4genes.png", 9, 7)
} else omitida("fig_KM_4genes.png", "falta resultados/cache_G_GSE62564.rds (script 04), GSE62564_series_matrix.txt o el paquete survival")

S("\nFiguras en: ", .LOG$dir)
end_log()
