# =============================================================================
# 10_diagrama_flujo.R   DIAGRAMA DE FLUJO DEL ANALISIS (Figura 1 de la memoria)
#   Dibuja con ggplot2 el flujo real de los scripts 01-09: datos -> preprocesado -> contrastes ->
#   interpretacion funcional -> programas y supervivencia -> validacion externa -> integracion y
#   priorizacion de genes. El texto de las cajas esta escrito aqui; si cambia el analisis, se cambia aqui.
# Salidas  : resultados/figuras/fig_diagrama_flujo.png
# Requiere : ggplot2
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need("ggplot2")
library(ggplot2)
start_log("diagrama")

# cajas: x, y = centro; w, h = tamano; fill = color; txt = texto (\n = salto)
box <- function(id, x, y, w, h, txt, fill = "#F2F2F2", face = "plain")
  data.frame(id = id, x = x, y = y, w = w, h = h, txt = txt, fill = fill, face = face)
COL_D <- "#DCE9F5"; COL_P <- "#EAEAEA"; COL_A <- "#E3F1E0"; COL_V <- "#FBE5D6"; COL_I <- "#FFF2CC"

B <- rbind(
  # fila 1: datos
  box("d1", 2.0, 11.2, 3.6, 1.5, "GSE18271\nmicroarray Affymetrix\n22 tumores primarios\n(12 GN, 10 GNB)", COL_D),
  box("d2", 6.0, 11.2, 3.6, 1.5, "GSE212877\nRNA-seq\n8 l\u00edneas celulares\nx 3 r\u00e9plicas", COL_D),
  box("d3", 10.0, 11.2, 3.6, 1.5, "GSE62564 (+ GSE49711)\nRNA-seq SEQC\n498 pacientes con\nsupervivencia", COL_D),
  box("d4", 14.0, 11.2, 3.6, 1.5, "TARGET-NBL\nRNA-seq (GDC)\n153 tumores primarios\n(validaci\u00f3n externa)", COL_V),
  # fila 2: preprocesado
  box("p1", 2.0, 9.0, 3.6, 1.5, "Control de calidad, RMA,\nsonda \u2192 gen (sonda de\nmayor expresi\u00f3n media):\n17.654 genes", COL_P),
  box("p2", 6.0, 9.0, 3.6, 1.5, "filterByExpr, TMM, voom\n15.627 genes", COL_P),
  box("p3", 10.0, 9.0, 3.6, 1.5, "Transcrito \u2192 gen (suma de\nRPM): 16.300 genes;\nMYCN y estadio de GSE49711", COL_P),
  box("p4", 14.0, 9.0, 3.6, 1.5, "Cuentas \u2192 TPM,\nfiltro de expresi\u00f3n:\n24.463 genes", COL_P),
  # fila 3: contrastes
  box("c1", 2.0, 6.8, 3.6, 1.5, "PCA y limma:\nGNB frente a GN", COL_A),
  box("c2", 6.0, 6.8, 3.6, 1.5, "PCA y limma:\nS frente a N pareado,\nl\u00ednea como efecto aleatorio", COL_A),
  box("c3", 10.0, 6.8, 3.6, 1.5, "limma: HR frente a non-HR,\nsin y con ajuste por MYCN", COL_A),
  # fila 4: interpretacion
  box("f1", 6.0, 4.6, 11.6, 1.3, "Enriquecimiento funcional (ORA GO y KEGG, GSEA) y composici\u00f3n celular (marcadores)", COL_A),
  # fila 5: programas y supervivencia
  box("s1", 10.0, 2.6, 3.6, 1.7, "Programas GO (CICLO, RIBO,\nSINAPSIS, INMUNE) y Cox\najustado por edad, MYCN\ny estadio en SEQC", COL_A),
  box("s2", 14.0, 2.6, 3.6, 1.7, "Mismos programas y mismo\nmodelo, sin reajustar,\nen TARGET-NBL", COL_V),
  # fila 6: integracion
  box("i1", 8.0, 0.5, 12.0, 1.3, "Integraci\u00f3n: concordancia entre datasets y selecci\u00f3n posterior de genes\n(descubrimiento en SEQC con FDR < 0,05; r\u00e9plica en TARGET-NBL)", COL_I)
)

# flechas (x0, y0, x1, y1)
A <- do.call(rbind, lapply(list(
  c(2, 10.45, 2, 9.75), c(6, 10.45, 6, 9.75), c(10, 10.45, 10, 9.75), c(14, 10.45, 14, 9.75),
  c(2, 8.25, 2, 7.55), c(6, 8.25, 6, 7.55), c(10, 8.25, 10, 7.55),
  c(2, 6.05, 4.2, 5.25), c(6, 6.05, 6, 5.25), c(10, 6.05, 7.8, 5.25),
  c(8.6, 3.95, 9.3, 3.45), c(14, 8.25, 14, 3.45), c(11.8, 2.6, 12.2, 2.6),
  c(10, 1.75, 9.2, 1.15), c(14, 1.75, 11.5, 1.15), c(5, 3.95, 6.2, 1.15)
), function(v) data.frame(x = v[1], y = v[2], xend = v[3], yend = v[4])))

p <- ggplot() +
  geom_rect(data = B, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill),
            colour = "grey35", linewidth = 0.4) +
  geom_text(data = B, aes(x = x, y = y, label = txt), size = 3.1, lineheight = 0.95) +
  geom_segment(data = A, aes(x = x, y = y, xend = xend, yend = yend), colour = "grey30", linewidth = 0.5,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  scale_fill_identity() + coord_equal(xlim = c(0, 16), ylim = c(-0.3, 12.1), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(4, 4, 4, 4))
dir.create(file.path(RES, "figuras"), showWarnings = FALSE)
ggsave(file.path(RES, "figuras", "fig_diagrama_flujo.png"), p, width = 9, height = 7.3, dpi = 200, bg = "white")
S("figura: fig_diagrama_flujo.png")
end_log()
