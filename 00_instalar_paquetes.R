# =============================================================================
# 00_instalar_paquetes.R  -  ejecutar UNA vez. Instala solo lo que falte.
# =============================================================================
cran <- c("ggplot2", "ggrepel", "pheatmap", "survival", "survminer", "data.table", "jsonlite")
bioc <- c("affy", "hgu133plus2.db", "hgu133plus2cdf", "limma", "edgeR", "AnnotationDbi", "org.Hs.eg.db",
          "clusterProfiler", "enrichplot", "DOSE")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
miss_c <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_c)) install.packages(miss_c)
miss_b <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_b)) BiocManager::install(miss_b, ask = FALSE, update = FALSE)
cat("\nPaquetes que siguen sin estar instalados: ",
    paste(c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)], collapse = ", "), "\n")
cat("\nVersiones:\n"); for (p in c(cran, bioc)) if (requireNamespace(p, quietly = TRUE)) cat(sprintf("  %-18s %s\n", p, as.character(packageVersion(p))))
