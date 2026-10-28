# TFM-extra: reanálisis crítico y multicohorte de la transcriptómica del neuroblastoma

Código en R del Trabajo Fin de Máster de **Arturo Fernández García** (Máster Universitario en Bioinformática, Universidad Europea; tutor: Dr. Giancarlo Castellano), convocatoria extraordinaria 2026.

El trabajo reanaliza de forma reproducible tres conjuntos públicos de GEO (GSE18271, GSE212877 y GSE62564/GSE49711) y valida los resultados en TARGET-NBL (GDC). Todas las cifras, tablas y figuras de la memoria salen de estos scripts.

## Declaración sobre el uso de IA
Parte del código se escribió y depuró con apoyo de la herramienta Claude (Anthropic), tal y como se declara en el apartado 8 de la memoria. El autor decidió el enfoque y los criterios de análisis, ejecutó todos los scripts sobre los datos originales y revisó los resultados, y asume la responsabilidad del contenido.

## Contenido
| Script | Función |
|---|---|
| `00_config.R` | Configuración común (rutas, umbrales, semilla, funciones de registro y enriquecimiento). Lo cargan los demás con `source()`. |
| `00_instalar_paquetes.R` | Instala los paquetes que falten (una sola vez). |
| `01_GSE18271_microarray.R` | Tumores primarios: RMA, control de calidad, colapso de sondas a genes, GNB frente a GN, marcadores de Schwann y neuronales, enriquecimiento. |
| `02_GSE212877_lineas_celulares.R` | Líneas celulares: filtrado, TMM, voom; contrastes A a F, sensibilidad del contraste D, enriquecimiento. |
| `03_GSE62564_cohorte_clinica.R` | Cohorte SEQC: agregación de transcritos a genes, HR frente a non-HR, Kaplan-Meier y Cox. |
| `04_GSE62564_ajuste_MYCN_estadio.R` | Ajuste por MYCN y estadio, modelos de Cox M0 a M4, genes concordantes entre conjuntos. |
| `05_puntuaciones_programas.R` | Programas biológicos definidos por Gene Ontology, Cox ajustado, riesgos proporcionales, terciles. |
| `06_validacion_TARGET_NBL.R` | Validación externa de los programas en TARGET-NBL, sin reajustar. |
| `07_figuras_finales.R` | Figuras finales de resultados y de concordancia. |
| `08_priorizacion_y_replica_genes.R` | Cribado de Cox de todos los genes comunes, replicación en TARGET-NBL, enriquecimiento y solapamientos. |
| `09_paneles_memoria.R`, `10_diagrama_flujo.R` | Paneles compuestos y diagrama de flujo de la memoria. |
| `99_ejecutar_todo.R` | Ejecuta 01 a 05 seguidos y junta los resúmenes (06 depende de datos externos). |
| `traduccion_GO.csv` | Traducción al castellano de los términos de GO que aparecen en las figuras. |

## Datos de entrada (no incluidos)
Los datos no se incluyen por tamaño y licencia. Hay que descargarlos y colocarlos en `Datos/`:
- GEO: GSE18271 (archivos CEL), GSE212877 (matriz de cuentas), GSE62564 y GSE49711 (series matrix).
- TARGET-NBL: Genomic Data Commons (Sample Sheet, datos clínicos y archivos `star_gene_counts.tsv`, acceso abierto).

## Cómo ejecutarlo
1. Define la carpeta de trabajo con la variable de entorno `TFM_ROOT` o edita `ROOT` al principio de cada script. Debe contener `scripts/`, `Datos/` y se creará `resultados/`.
2. Ejecuta `00_instalar_paquetes.R` una vez.
3. Ejecuta los scripts por orden: 01, 02, 03, 04, 05, 06, 07, 08, 09 y 10 (el 99 encadena del 01 al 05). El 08 necesita antes el 04 y el 06.
4. Cada análisis escribe en `resultados/` tablas CSV, figuras PNG, un `RESUMEN_*.txt` con las cifras principales, un registro y `sessionInfo.txt` con las versiones de los paquetes.

## Reproducibilidad
La semilla aleatoria está fijada (`SEED` en `00_config.R`). Umbrales: FDR (Benjamini-Hochberg) < 0,05 y |log2FC| > 1 para las listas de genes diferenciales. Los criterios de selección de genes del script 08 se fijaron antes de ejecutarlo sobre TARGET-NBL.
