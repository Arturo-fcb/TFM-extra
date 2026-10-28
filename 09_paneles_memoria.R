# =============================================================================
# 09_paneles_memoria.R   PANELES DE LA MEMORIA (todas las figuras compuestas salen de aqui)
#   Une, en columnas verticales rotuladas A, B, ..., las figuras individuales que ya dejo el script 07
#   (y los PCA/volcanos de los scripts 01-04) sin tocar su contenido. Asi cada figura de la memoria se
#   puede regenerar con un solo comando, sin editores de imagen.
# Entradas : resultados/figuras/*.png  y  resultados/GSE212877/PCA_8_lineas.png
# Salidas  : resultados/paneles/panel_*.png
# Requiere : png, grid (incluido en R)
# Ejecutar DESPUES de 07_figuras_finales.R
# =============================================================================
ROOT <- Sys.getenv("TFM_ROOT", "C:/ARTURO/Arturo universidad/Master/UE/TFM/extraordinaria")
source(file.path(ROOT, "scripts", "00_config.R"))
need(c("png"))
library(grid)
start_log("paneles")

FIG <- file.path(RES, "figuras")

# Une las imagenes `files` en una sola columna; la anchura comun es la de la mas ancha.
# Cada imagen se reescala a esa anchura (sin deformarla) y lleva una letra en la esquina.
panel_vertical <- function(files, nombre, letras = LETTERS, ancho_px = 1900) {
  falta <- files[!file.exists(files)]
  if (length(falta)) { S("  OMITIDO ", nombre, ": faltan ", paste(basename(falta), collapse = ", ")); return(invisible(NULL)) }
  imgs <- lapply(files, png::readPNG)
  alt  <- vapply(imgs, function(im) ancho_px * nrow(im) / ncol(im), numeric(1))
  H    <- ceiling(sum(alt)); marg <- 90; W <- ancho_px + marg      # margen izquierdo para la letra del panel
  png(OUT(nombre), width = W, height = H, res = 200, bg = "white")
  grid.newpage()
  y_top <- H
  for (i in seq_along(imgs)) {
    yv <- 1 - (H - y_top) / H
    pushViewport(viewport(x = marg / W, y = yv, width = ancho_px / W, height = alt[i] / H, just = c("left", "top")))
    grid.raster(imgs[[i]], width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
    popViewport()
    pushViewport(viewport(x = 0, y = yv, width = marg / W, height = alt[i] / H, just = c("left", "top")))
    grid.text(letras[i], x = unit(0.5, "npc"), y = unit(1, "npc") - unit(10, "pt"), just = c("centre", "top"),
              gp = gpar(fontface = "bold", fontsize = 22, col = "black"))
    popViewport()
    y_top <- y_top - alt[i]
  }
  dev.off()
  S("  panel: ", nombre, " (", length(files), " figuras)")
  invisible(OUT(nombre))
}

SH("Paneles de la memoria")
# Fig. 2: exploracion y expresion diferencial en tejido primario
panel_vertical(file.path(FIG, c("fig_PCA_GSE18271_22_muestras.png", "fig_volcano_GSE18271.png")), "panel_GSE18271_PCA_volcano.png")
# Fig. 3: enriquecimiento funcional GSE18271
panel_vertical(file.path(FIG, c("fig_GO_GSE18271_mas_en_GNB.png", "fig_GO_GSE18271_mas_en_GN.png")), "panel_GSE18271_GO.png")
# Fig. 4: lineas celulares (PCA de las 8 lineas y contraste D)
panel_vertical(c(file.path(RES, "GSE212877", "PCA_8_lineas.png"), file.path(FIG, "fig_volcano_GSE212877_D.png")), "panel_GSE212877_PCA_volcano.png")
# Fig. 5: enriquecimiento funcional GSE212877 (contraste D)
panel_vertical(file.path(FIG, c("fig_GO_GSE212877_mas_en_S.png", "fig_GO_GSE212877_mas_en_N.png")), "panel_GSE212877_GO.png")
# Fig. 6: HR frente a no-HR, sin ajustar y ajustado por MYCN
panel_vertical(file.path(FIG, c("fig_volcano_GSE62564_HR.png", "fig_volcano_GSE62564_HR_ajustado_MYCN.png")), "panel_GSE62564_volcanos.png")
# Anexo: GO de HR ajustado por MYCN
panel_vertical(file.path(FIG, c("fig_GO_GSE62564_HR_ajustado_mas_en_HR.png", "fig_GO_GSE62564_HR_ajustado_menos_en_HR.png")), "panel_GSE62564_GO_HRajustado.png")
# Fig. programas: terciles (SEQC y TARGET) y forest plot
panel_vertical(file.path(FIG, c("fig_supervivencia_terciles_4programas.png", "fig_forest_4programas_SEQC_vs_TARGET.png")), "panel_programas_SEQC_TARGET.png", ancho_px = 2200)
# Anexo: CICLO por estratos y terciles
panel_vertical(file.path(FIG, c("fig_forest_CICLO_estratos.png", "fig_supervivencia_terciles_CICLO.png")), "panel_CICLO_anexo.png")
end_log()
