# ==============================================================================
# PIPELINE Microarray GSE86282
# Hantavirus (HTNV) / HFRS — HUVEC
# Normalización ya hecha en GEO
# ==============================================================================

R.version.string
rm(list=ls())
getwd()

# ====================================
# BLOQUE 0: CONFIGURACIÓN Y LIBRERÍAS
# ====================================

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(stringr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Paths
raw_matrix_path <- "data/raw/human/GSE86282_RAW/GSE86282_series_matrix.txt.gz"
output_dir      <- "figures/human/GSE86282/"
results_dir     <- "results/human/GSE86282/"

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/human/", recursive = TRUE, showWarnings = FALSE)

# =========================================
# BLOQUE 1: CARGA DE DATOS Y METADATOS GEO
# =========================================

# Desde GEO directamente (error al intentar leer archivo descargado)
gse <- getGEO("GSE86282", GSEMatrix = TRUE)
eset <- gse[[1]]
data <- exprs(eset) #matriz limpia
metadata_clinica <- pData(eset)
head(data[, 1:3])


# Crear metadatos manualmente (simples, C vs V)
conditions <- c(rep("mock", 3), rep("infected", 3))

# Crear el coldata de forma directa
coldata <- data.frame(
  row.names = colnames(data),
  infection = factor(conditions, levels = c("mock", "infected"))
)
print(table(coldata$infection))

# ==============================================================================
# BLOQUE 2: QC + NORMALIZACIÓN + PCA
# ==============================================================================

# Los datos ya están en log2 y normalizados (quantile por Agilent)
# Verificar:
cat("Valores: mín=", min(data), ", máx=", max(data), "\n")

# Crear objeto ExpressionSet para limma
expr_matrix <- as.matrix(data)
eset <- ExpressionSet(assayData = expr_matrix,
                      phenoData = AnnotatedDataFrame(coldata))

# Filtrado de baja expresión
# Umbral: mantener si al menos 3 muestras (1 grupo completo) > median global + 1*IQR

expr_threshold <- median(as.matrix(eset)) + IQR(as.matrix(eset))
keep <- rowSums(exprs(eset) > expr_threshold) >= 3
eset_filtered <- eset[keep, ]

cat(sprintf("después de filtrado: %d (de %d originales)\n", 
            nrow(eset_filtered), nrow(eset)))

# PCA
expr_norm <- exprs(eset_filtered)
pca_data <- prcomp(t(expr_norm), scale. = TRUE)
pca_df <- as.data.frame(pca_data$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(coldata %>% rownames_to_column("sample"), by = "sample")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = infection, size = 3)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("mock" = "#74add1", "infected" = "#d73027")) +
  labs(title = "PCA: Mock vs Infected (GSE86282)",
       x = paste0("PC1 (", round(summary(pca_data)$importance[2,1]*100), "%)"),
       y = paste0("PC2 (", round(summary(pca_data)$importance[2,2]*100), "%)")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(paste0(output_dir, "PCA_plot.jpeg"),
       p_pca, width = 8, height = 6, quality = 95)
print(p_pca)

# ==============================================================================
# BLOQUE 3: MODELO LINEAL CON LIMMA
# ==============================================================================

# Crear matriz de diseño
design <- model.matrix(~ 0 + infection, data = coldata)
print(design)

# Ajuste del modelo lineal
fit <- lmFit(eset_filtered, design)

# Empirical Bayes (shrinkage de varianzas)
fit_eb <- eBayes(fit)

# Extraer resultados del contraste: infected - mock
contrast_matrix <- makeContrasts(infectioninfected - infectionmock, levels = design)
fit_contrast <- contrasts.fit(fit_eb, contrast_matrix)
fit_contrast <- eBayes(fit_contrast)

# Tabla de resultados
res <- topTable(fit_contrast, adjust.method = "BH", number = Inf) %>%
  rownames_to_column("probe_id")

cat(sprintf("Resultados: %d probes con estadísticas\n", nrow(res)))

# ==============================================================================
# BLOQUE 4: CLASIFICACIÓN DE DEGs Y VOLCANO PLOT
# ==============================================================================

#  Clasificar DEGs (logFC > 1, adj.P.Val < 0.05)
res_df <- res %>%
  mutate(
    DEG_class = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE                           ~ "NS"
    ),
    DEG_class = factor(DEG_class, levels = c("Up", "Down", "NS"))
  )

print(table(res_df$DEG_class))

# Top genes para etiquetas
top_labels <- res_df %>%
  filter(DEG_class != "NS") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

# 3. Volcano plot
volcano_plot <- ggplot(res_df, aes(x = logFC,
                                   y = -log10(adj.P.Val),
                                   color = DEG_class)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = "#d62728",
                                "Down" = "#1f77b4",
                                "NS"   = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(data = top_labels,
                           aes(label = probe_id),
                           size = 2.5, max.overlaps = 20,
                           color = "black") +
  labs(title = "Volcano plot — HTNV infected vs Mock (GSE86282)",
       x = "Log2 Fold Change",
       y = "-log10(adj.P.Val)",
       color = "DEG") +
  theme_minimal(base_size = 13)

ggsave(paste0(output_dir, "volcano_HTNV_vs_mock.jpeg"),
       volcano_plot, width = 8, height = 6, quality = 95)
print(volcano_plot)


# ==============================================================================
# BLOQUE 5: HEATMAP TOP 50 DEGs
# ==============================================================================

# Seleccionar top 50
top50_probes <- res_df %>%
  filter(DEG_class != "NS") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(probe_id)

# Matriz de expresión
mat_heatmap <- exprs(eset_filtered)[top50_probes, ]
mat_scaled <- t(scale(t(mat_heatmap)))

# Anotación de columnas
col_annotation <- coldata %>%
  dplyr::select(infection)

# Colores
ann_colors <- list(
  infection = c("mock" = "#74add1", "infected" = "#d73027")
)

# Heatmap
pheatmap(mat_scaled,
         annotation_col  = col_annotation,
         annotation_colors = ann_colors,
         show_rownames   = TRUE,
         show_colnames   = TRUE,
         fontsize_row    = 7,
         cluster_cols    = TRUE,
         cluster_rows    = TRUE,
         color = colorRampPalette(c("#313695", "white", "#a50026"))(100),
         main = "Top 50 DEGs — HTNV infected vs Mock (GSE86282)",
         filename = paste0(output_dir, "heatmap_top50_DEGs.jpeg"),
         width = 10, height = 12)

# ==============================================================================
# BLOQUE 6: ANÁLISIS FUNCIONAL (ORA: GO + KEGG)
# ==============================================================================
'''
Los identificadores de sonda de la plataforma GPL2182 (ASHGV...) no logran mapear a 
ENTREZID mediante las bases de datos de Bioconductor 
(org.Hs.eg.db, paquetes .db específicos o biomaRt). 
He agotado las vías estándar de mapeo sin obtener una correspondencia fiable para ejecutar 
enrichGO o enrichKEGG.
En ningún caso he sido capaz de obtener resultados
'''



#=================================================
# BLOQUE 7: VOLCANOS MIELOIDE vs LINFOIDE
#=================================================

# Coger gene sets de células inmunes de MSigDB C8 (cell type signatures)
immune_sets <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  filter(str_detect(gs_name, "MONOCYTE|NEUTROPHIL|MYELOID|BCELL|TCELL|NK_CELL|PLASMA"))

# Separar mieloides y linfoides automáticamente
genes_mieloides_auto <- immune_sets %>%
  filter(str_detect(gs_name, "MONOCYTE|NEUTROPHIL|MYELOID")) %>%
  pull(gene_symbol) %>% unique() %>%
  c("TLR3", "CGAS", "STING1", "FLT1", "CD55", "TLR8")  # Específico Hantavirus mieloides

genes_linfoides_auto <- immune_sets %>%
  filter(str_detect(gs_name, "BCELL|TCELL|NK_CELL|PLASMA")) %>%
  pull(gene_symbol) %>% unique() %>% 
  c("CDKN1A")  # Específico Hantavirus linfoide
# Cruzar con DEGs
res_mieloide <- res_df %>% filter(symbol %in% genes_mieloides_auto, !is.na(log2FoldChange))
res_linfoide <- res_df %>% filter(symbol %in% genes_linfoides_auto, !is.na(log2FoldChange))

cat("Genes mieloides:", nrow(res_mieloide), "\n")
cat("Genes linfoides:", nrow(res_linfoide), "\n")

# Volcano mieloide
plot_volcano(res_mieloide,
             "Volcano — Genes Mieloides (Acute vs Recovered)",
             paste0(output_dir, "volcano_mieloides_acute_vs_recovered.jpeg"))

# Volcano linfoide
plot_volcano(res_linfoide,
             "Volcano — Genes Linfoides / BCR (Acute vs Recovered)",
             paste0(output_dir, "volcano_linfoides_acute_vs_recovered.jpeg"))












# ==============================================================================
# BLOQUE 8: GUARDADO FINAL
# ==============================================================================


# Tabla completa de resultados (logFC, P.Value, adj.P.Val, etc.)
write.csv(res_df,
          file = paste0(results_dir, "limma_results_full.csv"),
          row.names = FALSE)

# Sondas significativamente expresadas (DEGs según tus filtros)
write.csv(res_df %>% filter(DEG_class != "NS"),
          file = paste0(results_dir, "limma_significant_probes.csv"),
          row.names = FALSE)

# Matriz de expresión procesada
norm_exprs <- exprs(eset_filtered)
write.csv(as.data.frame(norm_exprs),
          file = paste0(results_dir, "normalized_intensity_matrix.csv"))

