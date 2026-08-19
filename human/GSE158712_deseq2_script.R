# ==============================================================================
# PIPELINE bulk RNA-seq — GSE158712
# HFRS (Hemorrhagic Fever with Renal Syndrome) / Hantavirus
# Diseño: patient_group (severe/moderate) x time_period (acute/recovered)
# Herramienta principal: DESeq2 + GSEA
# Basado en: GSE_timepoints_deseq2_pipeline.R
#===============================================================================

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
library(GEOquery)
library(tibble)
library(dplyr)
library(stringr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ashr)
library(msigdbr)
library(patchwork)


geo_id          <- "GSE158712"
base_dir        <- "C:/Users/Jacobo/Desktop/BIOINFO/PRÁCTICAS/TFM"
raw_counts_path <- paste0(base_dir, "/data/raw/human/GSE158712_RAW/GSE158712_raw_counts_GRCh38.p13_NCBI.tsv")  # ⚠️ ajustar nombre exacto del TSV
output_dir      <- paste0(base_dir, "/figures/human/", geo_id, "/")
results_dir     <- paste0(base_dir, "/results/human/", geo_id, "/")

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(base_dir, "/data/processed/human/"), recursive = TRUE, showWarnings = FALSE)

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
counts_raw <- read.table(raw_counts_path, header = TRUE, row.names = 1, sep = "\t")

gse              <- getGEO(geo_id, GSEMatrix = TRUE)
metadata_clinica <- pData(gse[[2]])

# Las columnas relevantes son `patient group:ch1` y `time period:ch1`
coldata_tmp <- metadata_clinica %>%
  dplyr::select(
    geo_accession,
    patient_group = `patient group:ch1`,
    time_period   = `time period:ch1`
  ) %>%
  mutate(
    # Referencia: moderate (grupo menos severo) y acute (fase más activa)
    patient_group = factor(patient_group, levels = c("moderate", "severe")),
    time_period   = factor(time_period,   levels = c("acute", "recovered"))
  )

rownames(coldata_tmp) <- NULL
coldata <- coldata_tmp %>%
  column_to_rownames("geo_accession")

# Alineamiento muestras
muestras_comunes <- intersect(colnames(counts_raw), rownames(coldata))
data_final        <- counts_raw[, muestras_comunes]
coldata_final     <- coldata[muestras_comunes, , drop = FALSE]

# Eliminar muestras sin anotación
coldata_final <- coldata_final[!is.na(coldata_final$patient_group) & 
                                 !is.na(coldata_final$time_period), ]
data_final    <- data_final[, rownames(coldata_final)]

cat("Muestras alineadas:", ncol(data_final), "\n")
cat("\nDistribución del diseño:\n")
print(table(coldata_final$patient_group, coldata_final$time_period))


#======================================
# BLOQUE 2: QC + NORMALIZACIÓN (PCA)
#======================================

# diseño factorial — ambos factores en el modelo
dds <- DESeqDataSetFromMatrix(
  countData = round(data_final),
  colData   = coldata_final,
  design    = ~ patient_group + time_period + patient_group:time_period
)

# Filtrado de baja expresión
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# VST para PCA
vsd <- vst(dds, blind = FALSE)

# PCA coloreado por time_period, forma por patient_group
pca_data <- plotPCA(vsd, intgroup = c("time_period", "patient_group"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 color = time_period,
                                 shape = patient_group)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = c("acute" = "#d73027", "recovered" = "#4575b4"),
                     labels = c("acute" = "Aguda", "recovered" = "Recuperado")) +
  scale_shape_manual(values = c("moderate" = 16, "severe" = 17),
                     labels = c("moderate" = "Moderado", "severe" = "Severo")) +
  labs(title  = "PCA — GSE158712 HFRS",
       x      = paste0("PC1: ", pct_var[1], "% varianza"),
       y      = paste0("PC2: ", pct_var[2], "% varianza"),
       color  = "Fase clínica",
       shape  = "Severidad") +
  theme_minimal(base_size = 13)

save_pub(paste0(output_dir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = DPI_HIGH)
print(pca_plot)

# Ajuste del modelo
dds <- DESeq(dds)
cat("\nCoeficientes disponibles:\n")
print(resultsNames(dds))


#================================================================
# BLOQUE 3: CONTRASTES DIFERENCIALES
#  contrastes clínicamente relevantes para HFRS/BCR
#================================================================

# Función reutilizable para extraer + anotar resultados
get_results <- function(dds, contrast_vec, label) {
  res_s <- lfcShrink(dds,
                     contrast = contrast_vec,
                     type     = "ashr")
  res_s$symbol <- mapIds(org.Hs.eg.db,
                         keys      = rownames(res_s),
                         column    = "SYMBOL",
                         keytype   = "ENTREZID",
                         multiVals = "first")
  res_df <- as.data.frame(res_s) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(padj)) %>%
    mutate(
      DEG_class = case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE                               ~ "NS"
      ),
      DEG_class = factor(DEG_class, levels = c("Up", "Down", "NS")),
      contrast  = label
    )
  cat("\n→", label, "\n")
  print(table(res_df$DEG_class))
  return(res_df)
}

# Contraste 1: acute vs recovered (efecto fase clínica, independiente de severidad)
res_time <- get_results(dds,
                        contrast_vec = c("time_period", "acute", "recovered"),
                        label        = "acute_vs_recovered")

# Contraste 2: severe vs moderate (efecto severidad, independiente de fase)
res_sev  <- get_results(dds,
                        contrast_vec = c("patient_group", "severe", "moderate"),
                        label        = "severe_vs_moderate")
res_sev_enrich <- res_sev %>%
  mutate(DEG_class = case_when(
    padj < 0.1 & log2FoldChange >  0.5 ~ "Up",
    padj < 0.1 & log2FoldChange < -0.5 ~ "Down",
    TRUE                                ~ "NS"
  ),
  DEG_class = factor(DEG_class, levels = c("Up", "Down", "NS")))

cat("Severo relajado:\n")
print(table(res_sev_enrich$DEG_class))

# Guardar resultados
res_list <- list(
  acute_vs_recovered = res_time,
  severe_vs_moderate = res_sev
)

for (nm in names(res_list)) {
  write.csv(res_list[[nm]],
            paste0(results_dir, "DESeq2_", nm, "_full.csv"),    row.names = FALSE)
  write.csv(res_list[[nm]] %>% filter(DEG_class != "NS"),
            paste0(results_dir, "DESeq2_", nm, "_DEGs.csv"),    row.names = FALSE)
}


#===========================================
# BLOQUE 4: VOLCANO PLOTS (por contraste)
#===========================================


plot_volcano <- function(res_df, title_str, output_path) {
  
  lab_up   <- res_df %>% filter(DEG_class == "Up")   %>% arrange(padj) %>% slice_head(n = 15)
  lab_down <- res_df %>% filter(DEG_class == "Down") %>% arrange(padj) %>% slice_head(n = 15)
  top_labels <- bind_rows(lab_up, lab_down)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = DEG_class)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = ifelse(is.na(symbol), gene_id, symbol)),
      size = 2.5, max.overlaps = 30, color = "black",
      box.padding = 0.35, min.segment.length = 0
    ) +
    labs(title = title_str, x = "Log2 Fold Change", y = "-log10(padj)", color = "DEG") +
    theme_minimal(base_size = 13)
  
  return(p)
}

p_time <- plot_volcano(res_time, "Aguda vs Recuperado (HFRS)",
                       paste0(output_dir, "volcano_acute_vs_recovered.jpeg"))

p_sev  <- plot_volcano(res_sev,  "Severo vs Moderado (HFRS)",
                       paste0(output_dir, "volcano_severe_vs_moderate.jpeg"))

# Panel 2: Aguda vs Recuperado | Severo vs Moderado
p_panel_contrastes <- (p_time | p_sev) +
  plot_annotation(
    title    = "Contrastes principales — GSE158712 HFRS",
    subtitle = "Izquierda: efecto fase clínica  |  Derecha: efecto severidad",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )
save_pub(paste0(output_dir, "volcano_panel_contrastes.png"),
         p_panel_contrastes, width = 16, height = 6, dpi = DPI_STD)
print(p_panel_contrastes)


#========================================
# BLOQUE 5: HEATMAP TOP 50
#  anotación doble (severidad + fase)
#========================================

# Usar el contraste más informativo para BCR: acute vs recovered
top50 <- res_time %>%
  filter(DEG_class != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 50)

mat_heatmap <- assay(vsd)[top50$gene_id, ]
rownames(mat_heatmap) <- make.unique(as.character(top50$symbol))
mat_scaled  <- t(scale(t(mat_heatmap)))

# anotación con dos variables (recodificada a español solo para visualización;
# coldata_final mantiene los niveles originales para el diseño de DESeq2)
col_annotation <- coldata_final %>%
  dplyr::select(time_period, patient_group) %>%
  mutate(
    time_period   = recode(time_period,   "acute" = "Aguda", "recovered" = "Recuperado"),
    patient_group = recode(patient_group, "moderate" = "Moderado", "severe" = "Severo")
  )

ann_colors <- list(
  time_period   = c("Aguda"    = "#d73027", "Recuperado" = "#4575b4"),
  patient_group = c("Moderado" = "#74c476", "Severo"     = "#e6550d")
)

pheatmap(mat_scaled,
         annotation_col    = col_annotation,
         annotation_colors = ann_colors,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 7,
         cluster_cols      = TRUE,
         cluster_rows      = TRUE,
         color    = colorRampPalette(c("#313695", "white", "#a50026"))(100),
         main     = "Top 50 DEGs — Aguda vs Recuperado (HFRS)",
         filename = paste0(output_dir, "heatmap_top50_acute_vs_recovered.png"),
         width = 10, height = 12, res = DPI_STD)


#=================================================
# BLOQUE 6: ORA (GO BP + KEGG)
#=================================================

# Preparar listas de genes por contraste y dirección

gene_lists <- list(
  "Aguda Up"    = res_time      %>% filter(DEG_class == "Up")   %>% pull(gene_id),
  "Aguda Down"  = res_time      %>% filter(DEG_class == "Down") %>% pull(gene_id),
  "Severo Up"   = res_sev_enrich %>% filter(DEG_class == "Up")   %>% pull(gene_id),
  "Severo Down" = res_sev_enrich %>% filter(DEG_class == "Down") %>% pull(gene_id)
)
# Universo combinado
universe_all <- union(res_time$gene_id, res_sev$gene_id)

# GO BP comparativo con simplify (elimina redundancia semántica)
# GO BP con simplify — aplicar enrichGO+simplify por grupo, luego compareCluster
# Primero simplify individual, luego comparar
go_simp_list <- lapply(names(gene_lists), function(nm) {
  genes <- gene_lists[[nm]]
  if (length(genes) < 5) return(NULL)
  ego <- enrichGO(gene          = genes,
                  universe      = universe_all,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  ego_s <- clusterProfiler::simplify(ego, cutoff = 0.6, by = "p.adjust")
  df <- as.data.frame(ego_s)
  if (nrow(df) == 0) return(NULL)
  df$Cluster <- nm
  return(df)
})

go_simp_df <- bind_rows(go_simp_list) %>%
  mutate(
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]; as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    Cluster = factor(Cluster, levels = names(gene_lists))
  )

# Top 8 términos por cluster para no saturar
top_terms <- go_simp_df %>%
  group_by(Cluster) %>%
  slice_min(p.adjust, n = 8) %>%
  ungroup() %>%
  pull(Description) %>% unique()

go_plot_df <- go_simp_df %>% filter(Description %in% top_terms)

# GRÁFICO GO BP 

go_cat_colors <- c(
  "Antiviral / IFN"           = "#3B6D11",
  "Inmunidad adaptativa"      = "#185FA5",
  "Regulación inmune"         = "#534AB7",
  "Inflamación / bacteriana"  = "#BA7517",
  "Ciclo celular / DNA"       = "#993C1D",
  "Vascular / hemostasia"     = "#993556",
  "Migración / adhesión"      = "#0F6E56"
)

# Función de asignación semántica
assign_go_category <- function(term) {
  t <- tolower(term)
  dplyr::case_when(
    str_detect(t, "virus|viral|interferon|ifn|antiviral")                              ~ "Antiviral / IFN",
    str_detect(t, "leukocyte|lymphocyte|adaptive|mhc|antigen|cell killing|immunity")  ~ "Inmunidad adaptativa",
    str_detect(t, "regulation of immune|cytokine|ubiquit|chemokine")                   ~ "Regulación inmune",
    str_detect(t, "bacter|inflammat|humoral|peptidyl|nitrosyl|respiratory burst|neutrophil") ~ "Inflamación / bacteriana",
    str_detect(t, "cell cycle|mitot|meioti|chromosome|chromatin|dna repl|dna geomet|centrom|segregat") ~ "Ciclo celular / DNA",
    str_detect(t, "coagul|hemosta|platelet|wound|vascular|endothel|smooth muscle|angiogen") ~ "Vascular / hemostasia",
    str_detect(t, "migrat|adhesion|locomotion|chemotax|mesenchym")                     ~ "Migración / adhesión",
    TRUE ~ "Regulación inmune"
  )
}

# Construir data.frame con los 4 clusters, añadir signo y faceta
go_lollipop_df <- go_simp_df %>%
  filter(Description %in% top_terms) %>%
  mutate(
    direction = ifelse(str_detect(Cluster, "Up"), "Up", "Down"),
    contrast  = ifelse(str_detect(Cluster, "Aguda"), "Aguda\nvs Recuperado", "Severo\nvs Moderado"),
    GeneRatio_signed = ifelse(direction == "Up", GeneRatio_num, -GeneRatio_num),
    category  = assign_go_category(Description)
  )

# Normalizar transparencia dentro de cada cluster
go_lollipop_df <- go_lollipop_df %>%
  group_by(Cluster) %>%
  mutate(
    padj_norm = (p.adjust - min(p.adjust)) / (max(p.adjust) - min(p.adjust) + 1e-10),
    dot_alpha = 1 - 0.55 * padj_norm
  ) %>%
  ungroup()

# Ordenar términos: primero Down (izquierda), luego Up (derecha), por GeneRatio
go_lollipop_df <- go_lollipop_df %>%
  arrange(contrast, GeneRatio_signed) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Límite X simétrico por faceta
x_max <- max(abs(go_lollipop_df$GeneRatio_signed), na.rm = TRUE) * 1.15

p_go_comp <- ggplot(go_lollipop_df,
                    aes(x     = GeneRatio_signed,
                        y     = Description,
                        size  = Count,
                        fill  = category,
                        alpha = dot_alpha)) +
  # Segmento lollipop
  geom_segment(aes(x = 0, xend = GeneRatio_signed,
                   y = Description, yend = Description,
                   color = category),
               linewidth = 0.45, show.legend = FALSE) +
  # Punto
  geom_point(shape = 21, color = "white", stroke = 0.35) +
  # Línea central
  geom_vline(xintercept = 0, linewidth = 0.55, color = "grey35") +
  # Facetas por contraste
  facet_wrap(~ contrast, scales = "free_x", ncol = 2) +
  # Etiquetas de dirección dentro de cada panel
  geom_text(data = data.frame(
    contrast = c("Aguda\nvs Recuperado", "Aguda\nvs Recuperado",
                 "Severo\nvs Moderado",  "Severo\nvs Moderado"),
    label    = c("← Down", "Up →", "← Down", "Up →"),
    x        = c(-x_max * 0.82,  x_max * 0.82, -x_max * 0.82,  x_max * 0.82),
    y        = c(0.6, 0.6, 0.6, 0.6),
    hjust    = c(0, 1, 0, 1)
  ),
  aes(x = x, y = y, label = label, hjust = hjust),
  inherit.aes = FALSE,
  size = 2.8, color = "grey45", fontface = "italic") +
  # Escalas
  scale_fill_manual(values  = go_cat_colors, name = "Categoría funcional") +
  scale_color_manual(values = go_cat_colors) +
  scale_size_continuous(range = c(2.5, 9), name = "Count") +
  scale_alpha_identity() +
  scale_x_continuous(labels = function(x) round(abs(x), 2),
                     limits = c(-x_max, x_max)) +
  # Tema
  labs(
    title    = "GO Biological Process — Aguda vs Recuperado  |  Severo vs Moderado",
    subtitle = "Umbrales: Aguda (padj<0.05, |LFC|>1) | Severo (padj<0.1, |LFC|>0.5) | simplify cutoff=0.6\nTamaño = Count  |  Intensidad = p.adjust  |  Color = categoría funcional",
    x = "GeneRatio", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 8, color = "grey45"),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 8),
    legend.position    = "right",
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.major.x = element_line(color = "grey88"),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(1.5, "lines")
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, alpha = 0.9)),
    size = guide_legend(override.aes = list(fill = "grey60", alpha = 0.8))
  )

n_terms_total <- length(unique(go_lollipop_df$Description))
fig_height    <- max(9, n_terms_total * 0.38 + 3)

save_pub(paste0(output_dir, "GO_BP_comparativo.png"),
         p_go_comp, width = 14, height = fig_height, dpi = DPI_STD)
print(p_go_comp)
# Fin gráfico GO BP comparativo

#=================================================
# BLOQUE 7: GENES B/PLASMABLAST DE INTERÉS (BCR)
# Filtro 4-fold change (log2FC > 2) para IgG/IgM
#=================================================

# Genes marcadores de respuesta B humoral
genes_bcr_interest <- c(
  # Activación B / Plasmablasts
  "CD19", "CD20", "MS4A1", "CD38", "SDC1", "PRDM1", "IRF4", "XBP1",
  "AICDA", "MKI67",
  # Inmunoglobulinas (cadenas pesadas, isotipos relevantes para BCR)
  "IGHG1", "IGHG2", "IGHG3", "IGHG4",
  "IGHM", "IGHA1", "IGHA2",
  # Cadenas ligeras
  "IGKC", "IGLC1", "IGLC2",
  # Coestimulación / activación
  "CD27", "CD80", "CD86", "CXCR5", "BCL6",
  #interés hantavirus
  "CD55", "PCDH1", "TLR3", "CGAS", "TLR8", "FLT1", "STING1"
)

# Extraer estos genes del contraste acute_vs_recovered
bcr_res <- res_time %>%
  filter(symbol %in% genes_bcr_interest) %>%
  arrange(desc(abs(log2FoldChange)))

cat("\n=== GENES B/PLASMABLAST — Aguda vs Recuperado ===\n")
print(bcr_res %>% dplyr::select(symbol, log2FoldChange, padj, DEG_class))

# Aplicar criterio 4-fold change (log2FC > 2) para IgG/IgM — CRITERIO BCR
bcr_4fc <- bcr_res %>%
  filter(
    str_detect(symbol, "^IGH[GM]|^IGHG"),  # IgG e IgM
    abs(log2FoldChange) > 2                  # 4-fold change
  )

cat("\n=== BCR IgG/IgM con |LFC| > 2 (4-fold change) ===\n")
print(bcr_4fc %>% dplyr::select(symbol, log2FoldChange, padj, DEG_class))

write.csv(bcr_res,  paste0(results_dir, "BCR_genes_acute_vs_recovered.csv"),  row.names = FALSE)
write.csv(bcr_4fc,  paste0(results_dir, "BCR_IgGIgM_4fold_filtered.csv"),     row.names = FALSE)

# Heatmap específico genes BCR (si hay suficientes)
if (nrow(bcr_res) >= 5) {
  mat_bcr <- assay(vsd)[bcr_res$gene_id[bcr_res$gene_id %in% rownames(vsd)], ]
  rownames(mat_bcr) <- make.unique(as.character(
    bcr_res$symbol[bcr_res$gene_id %in% rownames(vsd)]
  ))
  mat_bcr_sc <- t(scale(t(mat_bcr)))
  
  pheatmap(mat_bcr_sc,
           annotation_col    = col_annotation,
           annotation_colors = ann_colors,
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 9,
           cluster_cols      = TRUE,
           cluster_rows      = TRUE,
           color    = colorRampPalette(c("#313695", "white", "#a50026"))(100),
           main     = "Genes B/Plasmablast — Aguda vs Recuperado",
           filename = paste0(output_dir, "heatmap_BCR_genes.png"),
           width = 10, height = 8, res = DPI_STD)
}

#=================================================
# BLOQUE 8: VOLCANOS MIELOIDE vs LINFOIDE
#=================================================

# Crear res_mieloide y res_linfoide
immune_sets <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  filter(str_detect(gs_name, "MONOCYTE|NEUTROPHIL|MYELOID|BCELL|TCELL|NK_CELL|PLASMA"))

genes_mieloides_auto <- immune_sets %>%
  filter(str_detect(gs_name, "MONOCYTE|NEUTROPHIL|MYELOID")) %>%
  pull(gene_symbol) %>% unique() %>%
  c("TLR3", "CGAS", "STING1", "FLT1", "CD55", "TLR8")

genes_linfoides_auto <- immune_sets %>%
  filter(str_detect(gs_name, "BCELL|TCELL|NK_CELL|PLASMA")) %>%
  pull(gene_symbol) %>% unique() %>%
  c("CDKN1A")

res_mieloide <- res_time %>% filter(symbol %in% genes_mieloides_auto, !is.na(log2FoldChange))
res_linfoide <- res_time %>% filter(symbol %in% genes_linfoides_auto, !is.na(log2FoldChange))

cat("Genes mieloides:", nrow(res_mieloide), "\n")
cat("Genes linfoides:", nrow(res_linfoide), "\n")

# Volcanos individuales
p_miel <- plot_volcano(res_mieloide,
                       "Genes Mieloides (Aguda vs Recuperado)",
                       paste0(output_dir, "volcano_mieloides_acute_vs_recovered.jpeg"))

p_linf <- plot_volcano(res_linfoide,
                       "Genes Linfoides / BCR (Aguda vs Recuperado)",
                       paste0(output_dir, "volcano_linfoides_acute_vs_recovered.jpeg"))

# Panel mieloides | linfoides
p_panel_inmune <- (p_miel | p_linf) +
  plot_annotation(
    title    = "Poblaciones inmunes — Aguda vs Recuperado (GSE158712)",
    subtitle = "Izquierda: genes mieloides  |  Derecha: genes linfoides/BCR",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )
save_pub(paste0(output_dir, "volcano_panel_miel_linf.png"),
         p_panel_inmune, width = 16, height = 6, dpi = DPI_STD)
print(p_panel_inmune)


#==========================================
# BLOQUE 9: EXPORTACIÓN FINAL
#==========================================

# Counts normalizados
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts),
          paste0(results_dir, "normalized_counts.csv"))
