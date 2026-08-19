# ==============================================================================
# PIPELINE bulk RNA-seq — GSE133634
# Hantavirus (HTNV) / HFRS
# Diseño:  infected vs control 
# Herramienta principal: DESeq2
# ==============================================================================

R.version.string
rm(list=ls())
getwd()

#========================================
# BLOQUE 0: CONFIGURACIÓN Y LIBRERÍAS
#========================================

library(DESeq2)
library(factoextra)
library(ggplot2)
library(pheatmap)
library(edgeR)
library(GEOquery)
library(tibble)
library(dplyr)
library(stringr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(patchwork)


# Paths
raw_counts_path <- "data/raw/human/GSE133634_RAW/GSE133634_raw_counts_GRCh38.p13_NCBI.tsv"
output_dir      <- "figures/human/GSE133634/"
results_dir     <- "results/human/GSE133634/"

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/human/", recursive = TRUE, showWarnings = FALSE)

# Parámetros de resolución 
DPI_STD  <- 300   # estándar para la mayoría de figuras
DPI_HIGH <- 600   # figuras con texto pequeño 

# Guardado de figuras 
save_pub <- function(path, plot, width, height, dpi = DPI_STD) {
  
  path <- sub("\\.jpe?g$", ".png", path, ignore.case = TRUE)
  
  png(filename = path,
      width    = width,
      height   = height,
      units    = "in",
      res      = dpi,
      type     = "cairo",
      bg       = "white")
  
  tryCatch(
    print(plot),
    error = function(e) {
      dev.off()
      warning(sprintf("Error renderizando %s: %s", basename(path), e$message))
      return(invisible(NULL))
    }
  )
  dev.off()
  cat(sprintf("  OK: %s  [%g x %g in | %d dpi | PNG]\n",
              basename(path), width, height, dpi))
}


#============================================
# BLOQUE 1: CARGA DE DATOS Y METADATOS GEO
#============================================

# Carga y metadata
data <- read.table(raw_counts_path, header=TRUE, row.names=1, sep="\t")

gse <- getGEO("GSE133634", GSEMatrix = TRUE)
metadata_clinica <- pData(gse[[1]])
coldata_tmp <- metadata_clinica %>%
  dplyr::select(geo_accession, infeccion_raw = `characteristics_ch1.2`) %>%
  mutate(
    infeccion = case_when(
      str_detect(infeccion_raw, "HTNV infected") ~ "infected",
      str_detect(infeccion_raw, "mock infected") ~ "mock"
    ),
    infeccion = factor(infeccion, levels = c("mock", "infected"))
  )

# Eliminar los nombres de fila previos
rownames(coldata_tmp) <- NULL

# Columna en nombres de fila
coldata <- coldata_tmp %>%
  column_to_rownames("geo_accession")

# Alineamiento
muestras_comunes <- intersect(colnames(data), rownames(coldata))
data_final <- data[, muestras_comunes]
coldata_final <- coldata[muestras_comunes, , drop = FALSE]
cat("Muestras alineadas:", ncol(data_final), "\n")

print(table(coldata_final$infeccion))

#======================================
# BLOQUE 2: QC + NORMALIZACIÓN (PCA)
#======================================

# Crear el objeto dds
dds <- DESeqDataSetFromMatrix(countData = round(data_final),
                              colData = coldata_final,
                              design = ~ infeccion)
# Filtrado de baja expresión (eliminar ruido)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep,]

# DESeq normalization
dds <- DESeq(dds)

# Extraer resultados (con shrinkage de LFC)
res_shrink <- lfcShrink(dds,
                        coef = "infeccion_infected_vs_mock",
                        type = "apeglm")
# Normalización VST para PCA
vsd <- vst(dds, blind = FALSE)
# Datos usados por plotPCA
pca_data <- plotPCA(vsd, intgroup = "infeccion", returnData = TRUE)

# Etiquetas (usa los nombres de las muestras)
pca_data$sample <- rownames(pca_data)

percentVar <- round(100 * attr(pca_data, "percentVar"))

pca <- ggplot(pca_data,
              aes(PC1, PC2, color = infeccion, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5,
                  max.overlaps = Inf,
                  box.padding = 0.4) +
  scale_color_manual(values = c("mock" = "#2166AC", "infected" = "#B2182B"),
                     labels = c("mock" = "Control", "infected" = "Infectado"),
                     name   = "Estado") +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA: Control vs Infectado") +
  theme_classic(base_size = 13)

save_pub(paste0(output_dir, "pca_plot.png"),
         pca,
         width = 8,
         height = 6,
         dpi = DPI_HIGH)
# Extraer resultados 
res <- results(dds)
summary(res)

#===========================================
# BLOQUE 3: FILTRADO DE DEGS + VOLCANOPLOT
#===========================================
# Mapeo de los genes
res_shrink$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_shrink),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Convertir a data.frame
res_df <- as.data.frame(res_shrink) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj))

write.table(res_df, 
            file = paste0(results_dir, "Resultados_DEGs_.tsv"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)



# Clasificar DEGs
res_df <- res_df %>%
  mutate(
    DEG_class = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                               ~ "NS"
    ),
    DEG_class = factor(DEG_class, levels = c("Up", "Down", "NS"))
  )

cat("\n✓ DEGs totales (padj<0.05, |LFC|>1):\n")
print(table(res_df$DEG_class))

# Top 20 para etiquetas
top_labels <- res_df %>%
  filter(DEG_class != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 20)


# Volcano plot
plot_volcano <- function(res_df, title_str, output_path) {
  
  top_up <- res_df %>%
    filter(DEG_class == "Up") %>%
    arrange(padj) %>%
    slice_head(n = 15)
  
  top_down <- res_df %>%
    filter(DEG_class == "Down") %>%
    arrange(padj) %>%
    slice_head(n = 10)           # etiqueta TODOS los Down (hasta 10)
  
  top_labels <- bind_rows(top_up, top_down)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = DEG_class)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1),     linetype = "dashed", color = "black") +
    ggrepel::geom_text_repel(data = top_labels,
                             aes(label = ifelse(is.na(symbol), gene_id, symbol)),
                             size = 2.5, max.overlaps = 20, color = "black") +
    labs(title = title_str, x = "Log2 Fold Change", y = "-log10(padj)", color = "DEG") +
    theme_minimal(base_size = 13)
  
  print(p)
}

#========================================
# BLOQUE 4: HEATMAP TOP 50
#========================================

# Top 50 DEGs
top50_genes <- res_df %>%
  filter(DEG_class != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 50)

# Matriz VST
mat_heatmap <- assay(vsd)[top50_genes$gene_id, ] #extraer top 50

rownames(mat_heatmap) <- make.unique(as.character(top50_genes$symbol))

mat_scaled <- t(scale(t(mat_heatmap))) #escalar

# Anotación de columnas
col_annotation <- coldata_final %>%
  dplyr::select(infeccion)

# Colores
ann_colors <- list(
  infeccion = c("mock" = "#74add1", "infected" = "#d73027")
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
         main = "Top 50 DEGs — Infectado vs Control",
         filename = paste0(output_dir, "heatmap_top50_DEGs.png"),
         width = 10, height = 12, res = DPI_STD)



#=================================================
# BLOQUE 5: ANÁLISIS FUNCIONAL (ORA: GO + KEGG) (GSEA)
#=================================================


#1. Separar up y down regulated
genes_up   <- res_df %>% filter(DEG_class == "Up")   %>% pull(gene_id)
genes_down <- res_df %>% filter(DEG_class == "Down") %>% pull(gene_id)
genes_all  <- res_df %>% filter(DEG_class != "NS")   %>% pull(gene_id)

# 2. Universo: todos los genes testeados
universe_genes <- res_df$gene_id

# 3. GO Biological Process — genes UP
ego_up <- enrichGO(gene          = genes_up,
                   universe      = universe_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",  
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
ego_up <- clusterProfiler::simplify(ego_up, cutoff = 0.6, by = "p.adjust") 

# 4. GO Biological Process — genes DOWN
ego_down <- enrichGO(gene          = genes_down,
                     universe      = universe_genes,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)
ego_down <- clusterProfiler::simplify(ego_down, cutoff = 0.6, by = "p.adjust")  


# 5. KEGG — todos los DEGs
ekegg <- enrichKEGG(gene          = genes_all,
                    universe      = universe_genes,
                    organism      = "hsa",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)

# 6. Dotplot KEGG (Figura 15 del TFM)
p_kegg  <- dotplot(ekegg,    showCategory = 15, title = "KEGG — All DEGs") +
  theme_minimal()

save_pub(paste0(output_dir, "KEGG_all.png"),   p_kegg,  width = 9, height = 8, dpi = DPI_STD)

print(p_kegg)


# 7 GRÁFICO COMBINADO BIDIRECCIONAL GO BP (Up + Down)

go_cat_colors <- c(
  "Antiviral / IFN"          = "#3B6D11",
  "Inmunidad adaptativa"     = "#185FA5",
  "Regulación inmune"        = "#534AB7",
  "Inflamación / bacteriana" = "#BA7517",
  "Diferenciación epidérmica"= "#993C1D",
  "Ubiquitinación / MHC"    = "#993556"
)

# Función para asignar categoría según el nombre del término GO
assign_go_category <- function(term) {
  term_low <- tolower(term)
  dplyr::case_when(
    str_detect(term_low, "virus|viral|interferon|ifn|antiviral") ~
      "Antiviral / IFN",
    str_detect(term_low, "leukocyte|lymphocyte|adaptive|mhc|antigen|cell killing") ~
      "Inmunidad adaptativa",
    str_detect(term_low, "regulation of immune|regulation of response|cytokine|ubiquit") ~
      "Regulación inmune",
    str_detect(term_low, "bacter|inflammat|humoral|peptidyl|nitrosyl|respiratory burst|leukocyte migrat|aggregat|modulation of process") ~
      "Inflamación / bacteriana",
    str_detect(term_low, "keratinocyte|epiderm|filament|hair cycle|keratinoc") ~
      "Diferenciación epidérmica",
    str_detect(term_low, "ubiquit|polyubiquit") ~
      "Ubiquitinación / MHC",
    TRUE ~ "Inflamación / bacteriana"   # categoría por defecto para términos no clasificados
  )
}

# Extraer tablas de resultados
df_up   <- as.data.frame(ego_up)   %>% mutate(direction = "Up")
df_down <- as.data.frame(ego_down) %>% mutate(direction = "Down")

# Tomar los top N términos de cada dirección
N_TERMS <- 15

df_up_top   <- df_up   %>% arrange(p.adjust) %>% slice_head(n = N_TERMS)
df_down_top <- df_down %>% arrange(p.adjust) %>% slice_head(n = N_TERMS)

# Construir data.frame combinado
df_combined <- bind_rows(df_up_top, df_down_top) %>%
  mutate(
    # GeneRatio como numérico
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- as.numeric(strsplit(x, "/")[[1]])
      parts[1] / parts[2]
    }),
    # Signo: Up positivo, Down negativo
    GeneRatio_signed = ifelse(direction == "Up", GeneRatio_num, -GeneRatio_num),
    # Categoría funcional
    category = assign_go_category(Description),
    # Orden visual: Down arriba, Up abajo (separados), dentro de cada grupo ordenar por GeneRatio
    sort_key = ifelse(direction == "Up", GeneRatio_num, -GeneRatio_num)
  ) %>%
  arrange(direction, sort_key) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Escalar transparencia por p.adjust (dentro de cada dirección para comparabilidad)
df_combined <- df_combined %>%
  group_by(direction) %>%
  mutate(
    padj_norm  = (p.adjust - min(p.adjust)) / (max(p.adjust) - min(p.adjust) + 1e-10),
    dot_alpha  = 1 - 0.55 * padj_norm   # más significativo → más opaco
  ) %>%
  ungroup()

# Construir el gráfico
p_go_combined <- ggplot(df_combined,
                        aes(x    = GeneRatio_signed,
                            y    = Description,
                            size = Count,
                            fill = category,
                            alpha = dot_alpha)) +
  # Línea de lollipop
  geom_segment(aes(x = 0, xend = GeneRatio_signed,
                   y = Description, yend = Description,
                   color = category),
               linewidth = 0.5, show.legend = FALSE) +
  # Punto
  geom_point(shape = 21, color = "white", stroke = 0.4) +
  # Línea vertical en 0
  geom_vline(xintercept = 0, linewidth = 0.6, color = "grey30") +
  # Anotación Up / Down
  annotate("text", x =  max(abs(df_combined$GeneRatio_signed)) * 0.9,
           y = 0.4, label = "Up-regulated →",
           hjust = 1, size = 3, color = "grey40", fontface = "italic") +
  annotate("text", x = -max(abs(df_combined$GeneRatio_signed)) * 0.9,
           y = 0.4, label = "← Down-regulated",
           hjust = 0, size = 3, color = "grey40", fontface = "italic") +
  # Escalas
  scale_fill_manual(values = go_cat_colors,  name = "Categoría funcional") +
  scale_color_manual(values = go_cat_colors) +
  scale_size_continuous(range = c(3, 10), name = "Count") +
  scale_alpha_identity() +
  scale_x_continuous(labels = function(x) round(abs(x), 2)) +
  # Tema
  labs(
    title    = "GO Biological Process — Up & Down regulated (Infectado vs Control)",
    subtitle = "Tamaño = Count  |  Intensidad = p.adjust  |  Color = categoría funcional",
    x        = "GeneRatio",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, color = "grey45"),
    axis.text.y     = element_text(size = 8),
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.major.x = element_line(color = "grey88"),
    panel.grid.minor   = element_blank()
  ) +
  guides(
    fill  = guide_legend(override.aes = list(size = 5, alpha = 0.9)),
    size  = guide_legend(override.aes = list(fill = "grey60", alpha = 0.8))
  )

# Guardar y mostrar
n_terms_total <- nrow(df_combined)
fig_height    <- max(8, n_terms_total * 0.38 + 2)

save_pub(paste0(output_dir, "GO_BP_combined_bidirectional.png"),
         p_go_combined,
         width   = 13,
         height  = fig_height,
         dpi     = DPI_STD)

print(p_go_combined)


# Ver top genes
cat("TOP 20 GENES UPREGULADOS:\n")
print(res_df %>% 
        filter(DEG_class == "Up") %>% 
        arrange(padj) %>% 
        slice_head(n=20) %>%
        dplyr::select(symbol, gene_id, log2FoldChange, padj))

cat("TOP 8 GENES DOWNREGULADOS:\n")
print(res_df %>% 
        filter(DEG_class == "Down") %>% 
        arrange(padj) %>%
        dplyr::select(symbol, gene_id, log2FoldChange, padj))

#=================================================
# BLOQUE 6: VOLCANOS MIELOIDE vs LINFOIDE
# Separar poblaciones
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

# Cruzar con  DEGs
res_mieloide <- res_df %>% filter(symbol %in% genes_mieloides_auto, !is.na(log2FoldChange))
res_linfoide <- res_df %>% filter(symbol %in% genes_linfoides_auto, !is.na(log2FoldChange))

cat("Genes mieloides:", nrow(res_mieloide), "\n")
cat("Genes linfoides:", nrow(res_linfoide), "\n")

#VOLCANO MIELOIDE
p_miel_133 <- plot_volcano(res_mieloide,
                           "Genes Mieloides — GSE133634",
                           paste0(output_dir, "volcano_mieloides.jpeg"))
#vOLCANO LINFOIDE
p_linf_133 <- plot_volcano(res_linfoide,
                           "Genes Linfoides / BCR — GSE133634",
                           paste0(output_dir, "volcano_linfoide.jpeg"))

# Panel combinado 
if (!is.null(p_miel_133) && !is.null(p_linf_133)) {
  p_combinado_133 <- p_miel_133 | p_linf_133 +
    plot_annotation(
      title    = "Mieloides vs Linfoides — Infectado vs Control (GSE133634)",
      subtitle = "Izquierda: genes mieloides  |  Derecha: genes linfoides/BCR",
      theme = theme(
        plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )
  save_pub(
    paste0(output_dir, "volcano_combinado_miel_linf.png"),
    p_combinado_133, width = 16, height = 7, dpi = DPI_STD
  )
  print(p_combinado_133)
}

#==========================================
# BLOQUE 7: GUARDADO FINAL
#==========================================

# Tabla completa
write.csv(res_df,
          file = paste0(results_dir, "DESeq2_results_full.csv"),
          row.names = FALSE)

# Solo DEGs
write.csv(res_df %>% filter(DEG_class != "NS"),
          file = paste0(results_dir, "DESeq2_DEGs_significant.csv"),
          row.names = FALSE)

# Enriquecimiento
if (nrow(as.data.frame(ego_up))   > 0)
  write.csv(as.data.frame(ego_up),   paste0(results_dir, "GO_BP_up.csv"),   row.names=FALSE)
if (nrow(as.data.frame(ego_down)) > 0)
  write.csv(as.data.frame(ego_down), paste0(results_dir, "GO_BP_down.csv"), row.names=FALSE)
if (nrow(as.data.frame(ekegg))    > 0)
  write.csv(as.data.frame(ekegg),    paste0(results_dir, "KEGG_all.csv"),   row.names=FALSE)

# Counts normalizados
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts),
          file = paste0(results_dir, "normalized_counts.csv"))

