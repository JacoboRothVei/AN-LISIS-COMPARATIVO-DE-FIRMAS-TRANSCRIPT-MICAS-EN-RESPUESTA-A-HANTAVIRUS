# ==============================================================================
# PIPELINE scRNA-seq — GSE161354
# Basado en: quadbio/scRNAseq_analysis_vignette (He & Treutlein, ETH Zürich)
# ==============================================================================

R.version.string
rm(list=ls())
getwd()

#=====================================
# BLOQUE 0: CONFIGURACIÓN
#=====================================

options(Seurat.object.assay.version = "v3")

library(Seurat)
library(GEOquery)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(SingleR)
library(celldex)
library(BiocParallel)
library(ggplot2)
library(ggrepel)
library(ragg)

# Rutas
raw_data_path <- "data/raw/human/GSE161354_RAW/"
output_dir    <- "figures/GSE161354/"
results_dir   <- "results/GSE161354/"

dir.create(output_dir,             recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir,            recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/human/",recursive = TRUE, showWarnings = FALSE)

all_plots <- list()

# PARÁMETROS DE RESOLUCIÓN 
DPI_STD   <- 300   # estándar para la mayoría de figuras
DPI_HIGH  <- 400   # figuras pequeñas 
QUAL_JPEG <- 95    # calidad JPEG (0-100)

save_pub <- function(path, plot, width, height, dpi = DPI_STD) {
  
  # Asegura extensión .png
  path <- sub("\\.jpe?g$", ".png", path, ignore.case = TRUE)
  
  png(
    filename = path,
    width    = width,
    height   = height,
    units    = "in",
    res      = dpi,
    type     = "cairo",
    bg       = "white"
  )
  
  tryCatch(
    {
      print(plot)
    },
    error = function(e) {
      dev.off()
      warning(sprintf("Error renderizando %s: %s", basename(path), e$message))
      return(invisible(NULL))
    }
  )
  
  dev.off()
}



#================================================
# BLOQUE 1 — CARGA DE DATOS Y METADATOS GEO
#================================================

gse <- getGEO("GSE161354", GSEMatrix = TRUE)
metadata_clinica <- pData(gse[[1]]) %>%
  dplyr::select(title, geo_accession, everything())

subcarpetas <- list.dirs(raw_data_path, full.names = TRUE, recursive = TRUE)
cat(sprintf("Subcarpetas detectadas: %d\n", length(subcarpetas)))

lista_seurat <- list()
for (carpeta in subcarpetas) {
  nombre <- basename(carpeta)
  tryCatch({
    counts <- Read10X(data.dir = paste0(carpeta, "/"))
    lista_seurat[[nombre]] <- CreateSeuratObject(
      counts = counts, project = nombre, min.cells = 0, min.features = 200)
    cat(sprintf("  Cargado: %s\n", nombre))
  }, error = function(e) cat(sprintf("  Error %s: %s\n", nombre, e$message)))
}

if (length(lista_seurat) == 0) stop("No se cargó ningún dataset.")

seurat_obj <- if (length(lista_seurat) == 1) {
  lista_seurat[[1]]
} else {
  merge(lista_seurat[[1]], y = lista_seurat[-1], add.cell.ids = names(lista_seurat))
}
Project(seurat_obj) <- "Hanta_GSE161354"
seurat_obj$orig.ident <- sapply(strsplit(colnames(seurat_obj), "_"), `[`, 1)

metadata_geo_mini <- metadata_clinica %>%
  dplyr::select(title,
                edad           = `age:ch1`,
                estado_clinico = `disease state:ch1`,
                sexo           = `Sex:ch1`)

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(metadata_geo_mini, by = c("orig.ident" = "title")) %>%
  column_to_rownames("cell_id")

cat(sprintf("%d células | %d muestras\n",
            ncol(seurat_obj), length(unique(seurat_obj$orig.ident))))

#=====================================
# BLOQUE 2: CONTROL DE CALIDAD
#=====================================

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT[-\\.]")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

qc_por_id <- seurat_obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(n_celulas      = n(),
            mediana_genes  = round(median(nFeature_RNA)),
            mediana_counts = round(median(nCount_RNA)),
            mediana_mt     = round(median(percent.mt), 2),
            mediana_rb     = round(median(percent.rb), 2),
            pct_excluir_mt    = round(100 * mean(percent.mt > 15), 1),
            pct_excluir_genes = round(100 * mean(nFeature_RNA < 500 | nFeature_RNA > 4000), 1),
            .groups = "drop") %>%
  arrange(desc(mediana_mt))

write.csv(qc_por_id, paste0(results_dir, "QC_por_paciente.csv"), row.names = FALSE)

  
p_vln_pre <- VlnPlot(seurat_obj,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
                   ncol = 4, group.by = "orig.ident", pt.size = 0) +
plot_annotation(title = "QC pre-filtro — distribución por muestra")

p_scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                           feature2 = "percent.mt",   group.by = "orig.ident", pt.size = 0.3) +
ggtitle("% mitocondrial")
p_scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                           feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.3) +
ggtitle("Genes por célula")

p_scatter  <- p_scatter1 + p_scatter2 +
plot_annotation(title = "Scatter QC — diagnóstico de dobletes")

all_plots[["02_QC_violin"]]  <- p_vln_pre
all_plots[["03_QC_scatter"]] <- p_scatter

seurat_filtered <- subset(seurat_obj,
                        subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &
                          percent.mt < 15 & nCount_RNA < 25000)

cat(sprintf("Células iniciales: %d → retenidas: %d (%.1f%%)\n",
          ncol(seurat_obj), ncol(seurat_filtered),
          100 * ncol(seurat_filtered) / ncol(seurat_obj)))

saveRDS(seurat_filtered, "results/seurat_GSE161354_filtered.rds")


#========================================================
# BLOQUE 3: NORMALIZACIÓN
#========================================================

seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered, nfeatures = 3000)
seurat_filtered <- ScaleData(seurat_filtered)

#===================================================
# BLOQUE 4: REDUCCIÓN DIMENSIONAL Y CLUSTERING
#===================================================
seurat_filtered <- RunPCA(seurat_filtered, npcs = 50, verbose = FALSE)

#----------------------------------------------------------------
# PCA PSEUDOBULK POR PACIENTE
# Equivalente a plotPCA de DESeq2 — cada punto = 1 paciente
# Método: AggregateExpression → CPM log1p → prcomp top 2000 genes
#----------------------------------------------------------------
pseudobulk_counts <- AggregateExpression(
  seurat_filtered, group.by = "orig.ident",
  assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA

pseudobulk_counts <- pseudobulk_counts[rowSums(pseudobulk_counts) >= 10, ]

pseudobulk_log <- log1p(
  sweep(pseudobulk_counts, 2, colSums(pseudobulk_counts), "/") * 1e6)

vars_genes <- apply(pseudobulk_log, 1, var)
top_var    <- names(sort(vars_genes, decreasing = TRUE))[1:2000]
pca_res    <- prcomp(t(pseudobulk_log[top_var, ]), scale. = TRUE)
pct_var    <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

df_pca_pb <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("orig.ident") %>%
  left_join(seurat_filtered@meta.data %>%
              dplyr::select(orig.ident, estado_clinico) %>% distinct(),
            by = "orig.ident") %>%
  mutate(estado_clinico = factor(estado_clinico,
                                 levels = c("normal health control", "HFRS fever stage")))

p_pca_pb <- ggplot(df_pca_pb,
                   aes(x = PC1, y = PC2, color = estado_clinico, label = orig.ident)) +
  geom_point(size = 5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.5, color = "black",
                           box.padding = 0.4, max.overlaps = 20) +
  scale_color_manual(values = c("normal health control" = "#2166AC",
                                "HFRS fever stage"      = "#B2182B"),
                     name = "Estado clínico") +
  labs(title    = "PCA — GSE161354",
       subtitle = sprintf("Cada punto = 1 paciente  |  n=6 HFRS + n=2 Control\nPC1: %s%%  |  PC2: %s%%",
                          pct_var[1], pct_var[2]),
       x = sprintf("PC1: %s%% varianza", pct_var[1]),
       y = sprintf("PC2: %s%% varianza", pct_var[2])) +
  theme_minimal(base_size = 13) +
  theme(plot.subtitle = element_text(color = "grey40", size = 9),
        legend.position = "bottom")

all_plots[["04b_PCA_pseudobulk"]] <- p_pca_pb
cat("PCA pseudobulk generado.\n")

n_dims <- 15
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:n_dims)
seurat_filtered <- FindClusters(seurat_filtered,  resolution = 0.5)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:n_dims)

p_umap_clusters <- DimPlot(seurat_filtered, reduction = "umap",
                           label = TRUE, pt.size = 0.4) +
  ggtitle("UMAP — clusters Seurat (pre-anotación)")
all_plots[["05_UMAP_clusters"]] <- p_umap_clusters

saveRDS(seurat_filtered, "data/processed/human/seurat_GSE161354_processed.rds")


#===========================================
# BLOQUE 5: ANOTACIÓN CELULAR CON SingleR
#===========================================
if (!"celltype_singler" %in% colnames(seurat_filtered@meta.data)) {
  
  seurat_filtered <- NormalizeData(seurat_filtered, assay = "RNA")
  ref   <- celldex::HumanPrimaryCellAtlasData()
  datos <- GetAssayData(seurat_filtered, assay = "RNA", layer = "data")
  
  genes_comunes <- length(intersect(rownames(datos), rownames(ref)))
  cat(sprintf("   Genes en común con referencia: %d\n", genes_comunes))
  if (genes_comunes < 1000) stop("< 1000 genes comunes — revisar formato.")
  
  options(future.globals.maxSize = 3000 * 1024^2)
  pred <- SingleR(test = datos, ref = ref, labels = ref$label.main,
                  BPPARAM = SerialParam())
  
  seurat_filtered$celltype_singler <- pred$labels
  seurat_filtered$celltype_pruned  <- pred$pruned.labels
  print(table(seurat_filtered$celltype_singler))
  
  saveRDS(seurat_filtered, "data/processed/human/seurat_GSE161354_processed.rds")
}

seurat_filtered <- RunPCA(seurat_filtered, npcs = 30, verbose = FALSE)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:15)

p_umap_singler <- DimPlot(seurat_filtered, reduction = "umap",
                          group.by = "celltype_singler",
                          label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle("UMAP — Anotación SingleR (HumanPrimaryCell)")
all_plots[["06_UMAP_SingleR"]] <- p_umap_singler

#===============================================
# BLOQUE 6: DEGs EN CÉLULAS INMUNES MADURAS
#===============================================
progenitores  <- c("HSC_-G-CSF", "CMP", "GMP", "MEP", "Pro-Myelocyte",
                   "Pre-B_cell_CD34-", "Pro-B_cell_CD34+", "BM", "BM & Prog.")
seurat_inmune <- subset(seurat_filtered,
                        subset = celltype_singler %in% progenitores, invert = TRUE)

cat(sprintf("Células inmunes maduras: %d | Tipos: %s\n",
            ncol(seurat_inmune),
            paste(unique(seurat_inmune$celltype_singler), collapse = ", ")))

Idents(seurat_inmune)      <- "celltype_singler"
DefaultAssay(seurat_inmune) <- "RNA"

markers <- FindAllMarkers(seurat_inmune, only.pos = FALSE,
                          min.pct = 0.25, logfc.threshold = 0.25,
                          assay = "RNA", test.use = "wilcox")

write.csv(markers, paste0(results_dir, "DEGs_celulas_inmunes.csv"), row.names = FALSE)

# Proporciones celulares
seurat_inmune$estado_clinico <- factor(seurat_inmune$estado_clinico)
prop_table <- prop.table(table(seurat_inmune$celltype_singler,
                               seurat_inmune$estado_clinico), margin = 2)
df_prop    <- as.data.frame(prop_table)

p_prop <- ggplot(df_prop, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Proporción de tipos celulares por estado clínico",
       x = "Estado", y = "Frecuencia")

write.csv(df_prop, paste0(results_dir, "Proporciones_Celulares_Clinicas.csv"), row.names = FALSE)
all_plots[["07_Proporciones_clinicas"]] <- p_prop

#=====================================
# BLOQUE 7: TRIPLE FILTRO GO
#=====================================
genes_deg     <- unique(markers$gene)
genes_validos <- genes_deg[genes_deg %in% keys(org.Hs.eg.db, keytype = "SYMBOL")]
genes_entrez  <- bitr(genes_validos, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_mf <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "MF", readable = TRUE,
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

terminos_receptor <- paste(c("receptor", "binding", "signaling", "pattern recognition",
                             "integrin", "adhesion", "endocytosis", "viral", "defense",
                             "immune", "innate", "phagocytosis", "antigen"), collapse = "|")

genes_go_mf <- unique(unlist(strsplit(
  go_mf@result %>%
    dplyr::filter(str_detect(Description, regex(terminos_receptor, ignore_case = TRUE))) %>%
    pull(geneID), "/")))

go_cc <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "CC", readable = TRUE,
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

terminos_membrana <- paste(c("plasma membrane", "cell surface", "extracellular",
                             "membrane receptor", "external side",
                             "integral component of membrane"), collapse = "|")

genes_go_cc <- unique(unlist(strsplit(
  go_cc@result %>%
    dplyr::filter(str_detect(Description, regex(terminos_membrana, ignore_case = TRUE))) %>%
    pull(geneID), "/")))

genes_cpdb <- tryCatch({
  read_csv(paste0("https://raw.githubusercontent.com/",
                  "Teichlab/cellphonedb-data/master/data/gene_input.csv"),
           show_col_types = FALSE) %>% pull(gene_name)
}, error = function(e) { warning("CellPhoneDB no disponible"); genes_go_cc })

cat(sprintf("   GO-MF: %d | GO-CC: %d | CellPhoneDB: %d\n",
            length(genes_go_mf), length(genes_go_cc), length(genes_cpdb)))

genes_candidatos_simbolos <- bitr(
  Reduce(intersect, list(
    genes_entrez$ENTREZID[genes_entrez$SYMBOL %in% genes_go_mf],
    genes_entrez$ENTREZID[genes_entrez$SYMBOL %in% genes_go_cc],
    genes_entrez$ENTREZID[genes_entrez$SYMBOL %in% genes_cpdb])),
  fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL

cat(sprintf("Candidatos tras filtro: %d\n", length(genes_candidatos_simbolos)))

#=====================================
# BLOQUE 8: RANKING Y VALIDACIÓN
#=====================================
candidatos <- markers %>%
  dplyr::filter(gene %in% genes_candidatos_simbolos) %>%
  dplyr::mutate(score = avg_log2FC * (-log10(p_val_adj + 1e-300))) %>%
  arrange(desc(score))

if (nrow(candidatos) == 0) {
  warning("Triple filtro vacío — usando doble filtro GO-MF + GO-CC")
  candidatos <- markers %>%
    dplyr::filter(gene %in% intersect(genes_go_mf, genes_go_cc)) %>%
    dplyr::mutate(score = avg_log2FC * (-log10(p_val_adj + 1e-300))) %>%
    arrange(desc(score))
}

receptores_dirigidos <- c("CD55", "PCDH1", "CGAS", "TLR3", "TLR8", "TMEM173", "FLT1")
candidatos           <- candidatos %>% dplyr::mutate(es_dirigido = gene %in% receptores_dirigidos)

write.csv(candidatos, paste0(results_dir, "candidatos_receptores_final.csv"), row.names = FALSE)

top_genes        <- candidatos %>% distinct(gene, .keep_all = TRUE) %>%
  slice_max(score, n = 15) %>% pull(gene)
top_genes_val    <- top_genes[top_genes %in% rownames(seurat_inmune)]
recept_presentes <- receptores_dirigidos[receptores_dirigidos %in% rownames(seurat_filtered)]

#=====================================
# BLOQUE 9: VISUALIZACIÓN
#=====================================
p_umap_inmune <- DimPlot(seurat_inmune, reduction = "umap",
                         group.by = "celltype_singler",
                         label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle("UMAP — Células inmunes maduras (PBMCs HFRS)")
all_plots[["10_UMAP_inmune"]] <- p_umap_inmune

if (length(top_genes_val) > 0) {
  all_plots[["11_Feature_candidatos"]] <- FeaturePlot(
    seurat_inmune, features = top_genes_val,
    cols = c("lightgrey", "darkred"), ncol = 3,
    order = TRUE, min.cutoff = "q10", pt.size = 0.3) &
    theme(plot.title = element_text(size = 9))
}

if (length(top_genes_val) > 0) {
  all_plots[["12_Dot_candidatos"]] <- DotPlot(
    seurat_inmune, features = top_genes_val,
    group.by = "celltype_singler", cols = c("lightgrey", "darkred")) +
    RotatedAxis() + ggtitle("Candidatos receptores por tipo celular — no sesgado")
}

if (length(recept_presentes) > 0) {
  all_plots[["13_Feature_dirigidos"]] <- FeaturePlot(
    seurat_filtered, features = recept_presentes,
    cols = c("lightgrey", "darkblue"), ncol = 2,
    order = TRUE, min.cutoff = 0, pt.size = 0.3) &
    theme(plot.title = element_text(size = 10, face = "bold"))
  
  recept_en_inmune <- recept_presentes[recept_presentes %in% rownames(seurat_inmune)]
  if (length(recept_en_inmune) > 0) {
    all_plots[["14_Dot_dirigidos"]] <- DotPlot(
      seurat_inmune, features = recept_en_inmune,
      group.by = "celltype_singler", cols = c("lightgrey", "darkblue")) +
      RotatedAxis() + ggtitle("Receptores dirigidos (lista tutora) por tipo celular")
  }
}

#=================================================
# BLOQUE 10: VOLCANOS MIELOIDE vs LINFOIDE
#=================================================
genes_hanta     <- c("TLR3", "CGAS", "STING1", "FLT1", "CD55", "TLR8", "CDKN1A")
tipos_mieloides <- c("Monocyte", "DC", "Macrophage", "Neutrophils")
tipos_linfoides <- c("T_cells", "NK_cell", "B_cell")

tipos_celulares <- unique(seurat_inmune$celltype_singler)
deg_por_tipo    <- list()

for (tipo in tipos_celulares) {
  subset_tipo <- subset(seurat_inmune, subset = celltype_singler == tipo)
  n_hfrs    <- sum(subset_tipo$estado_clinico == "HFRS fever stage",     na.rm = TRUE)
  n_control <- sum(subset_tipo$estado_clinico == "normal health control", na.rm = TRUE)
  if (n_hfrs < 10 || n_control < 10) {
    cat(sprintf("  %s: omitido\n", tipo)); next
  }
  Idents(subset_tipo) <- "estado_clinico"
  tryCatch({
    deg <- FindMarkers(subset_tipo,
                       ident.1 = "HFRS fever stage", ident.2 = "normal health control",
                       test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1,
                       assay = "RNA") %>%
      rownames_to_column("gene") %>%
      mutate(celltype  = tipo,
             DEG_class = factor(case_when(
               p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Up",
               p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
               TRUE ~ "NS"), levels = c("Up", "Down", "NS")),
             es_hanta  = gene %in% genes_hanta)
    deg_por_tipo[[tipo]] <- deg
  }, error = function(e) cat(sprintf("  Error %s: %s\n", tipo, e$message)))
}

deg_infeccion_all <- bind_rows(deg_por_tipo)
write.csv(deg_infeccion_all,
          paste0(results_dir, "DEGs_HFRS_vs_Control_por_celltipo.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# Función volcano 
# --------------------------------------------------------------------------
plot_volcano_hq <- function(deg_df, tipos, titulo) {
  
  df <- deg_df %>%
    filter(celltype %in% tipos) %>%
    group_by(gene) %>%
    slice_min(p_val_adj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(logp = -log10(p_val_adj + 1e-300))
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  # Top 15 UP por LFC 
  top_up <- df %>%
    filter(DEG_class == "Up") %>%
    slice_max(avg_log2FC, n = 15)
  
  # Top 15 DOWN por LFC (más negativos)
  top_down <- df %>%
    filter(DEG_class == "Down") %>%
    slice_min(avg_log2FC, n = 15)
  
  # Outliers extremos en Y
  top_y <- df %>% slice_max(logp, n = 5)
  
  # Genes hanta siempre etiquetados
  hanta_lab <- df %>% filter(es_hanta)
  
  top_lab <- bind_rows(top_up, top_down, top_y, hanta_lab) %>%
    distinct(gene, .keep_all = TRUE)
  
  ggplot(df, aes(x = avg_log2FC, y = logp, color = DEG_class)) +
    geom_point(data = df %>% filter(DEG_class == "NS"),
               alpha = 0.30, size = 1.6, shape = 16) +
    geom_point(data = df %>% filter(DEG_class != "NS"),
               alpha = 0.60, size = 1.8, shape = 16) +
    scale_color_manual(
      values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey72"),
      labels = c("Up" = "Up", "Down" = "Down", "NS" = "NS"),
      drop   = FALSE
    ) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed", color = "black", linewidth = 0.5) +
    ggrepel::geom_text_repel(
      data          = top_lab,
      aes(label     = gene),
      size          = 3.4,           
      color         = "black",
      fontface      = "bold",
      box.padding   = 0.45,
      point.padding = 0.3,
      max.overlaps  = 50,
      min.segment.length = 0.2,
      segment.size  = 0.3,
      segment.color = "grey50"
    ) +
    labs(
      title    = titulo,
      subtitle = "HFRS fever stage vs normal health control",
      x        = "Log2 Fold Change (HFRS / Control)",
      y        = expression(-log[10](adj.P.Val)),
      color    = "DEG"
    ) +
    theme_minimal(base_size = 14) +    # era 12
    theme(
      legend.position      = "bottom",
      legend.key.size      = unit(0.9, "lines"),
      legend.text          = element_text(size = 12),
      legend.title         = element_text(size = 12),
      plot.subtitle        = element_text(color = "grey40", size = 10),
      plot.title           = element_text(size = 14, face = "bold"),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(linewidth = 0.3, color = "grey88"),
      axis.title           = element_text(size = 13),
      axis.text            = element_text(size = 11)
    )
}

p_miel_hq <- plot_volcano_hq(deg_infeccion_all, tipos_mieloides, "Mieloides — HFRS vs Control")
p_linf_hq <- plot_volcano_hq(deg_infeccion_all, tipos_linfoides, "Linfoides — HFRS vs Control")

if (!is.null(p_miel_hq)) {
  ragg::agg_png(
    filename   = paste0(output_dir, "Fig12a_volcano_mieloide_HQ.png"),
    width      = 13, height = 9.5, units = "in",
    res        = DPI_HIGH, background = "white"
  )
  print(p_miel_hq)
  dev.off()
}

if (!is.null(p_linf_hq)) {
  ragg::agg_png(
    filename   = paste0(output_dir, "Fig12b_volcano_linfoide_HQ.png"),
    width      = 13, height = 9.5, units = "in",
    res        = DPI_HIGH, background = "white"
  )
  print(p_linf_hq)
  dev.off()
}

all_plots[["12a_volcano_mieloide"]] <- p_miel_hq
all_plots[["12b_volcano_linfoide"]] <- p_linf_hq

#================================================================
# BLOQUE 12: EXPRESIÓN DIFERENCIAL DE RECEPTORES CANDIDATOS
#================================================================
receptores_todos     <- unique(c(receptores_dirigidos, top_genes_val))
receptores_presentes <- receptores_todos[receptores_todos %in% rownames(seurat_inmune)]

deg_receptores_lista <- list()
for (tipo in unique(seurat_inmune$celltype_singler)) {
  subset_tipo <- subset(seurat_inmune, subset = celltype_singler == tipo)
  n_hfrs    <- sum(subset_tipo$estado_clinico == "HFRS fever stage",     na.rm = TRUE)
  n_control <- sum(subset_tipo$estado_clinico == "normal health control", na.rm = TRUE)
  if (n_hfrs < 10 || n_control < 10) { cat(sprintf("   %s: omitido\n", tipo)); next }
  Idents(subset_tipo) <- "estado_clinico"
  tryCatch({
    deg <- FindMarkers(subset_tipo,
                       ident.1 = "HFRS fever stage", ident.2 = "normal health control",
                       features = receptores_presentes, test.use = "wilcox",
                       min.pct = 0.05, logfc.threshold = 0, assay = "RNA") %>%
      rownames_to_column("gene") %>%
      mutate(celltype  = tipo,
             DEG_class = case_when(
               p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Up",
               p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
               TRUE ~ "NS"))
    deg_receptores_lista[[tipo]] <- deg
  }, error = function(e) cat(sprintf("   %s: error — %s\n", tipo, e$message)))
}

deg_receptores_all <- bind_rows(deg_receptores_lista)
write.csv(deg_receptores_all,
          paste0(results_dir, "DEG_receptores_por_tipo_celular.csv"), row.names = FALSE)

df_heat <- deg_receptores_all %>%
  mutate(sig_label        = case_when(p_val_adj < 0.001 ~ "***", p_val_adj < 0.01 ~ "**",
                                      p_val_adj < 0.05  ~ "*",   TRUE ~ ""),
         avg_log2FC_clamp = pmax(pmin(avg_log2FC, 3), -3))

p_heat_receptores <- ggplot(df_heat, aes(x = celltype, y = gene, fill = avg_log2FC_clamp)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.5, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-3, 3),
                       name = "LFC\n(HFRS/ctrl)\n[±3]", na.value = "grey85") +
  labs(title    = "Expresión diferencial de receptores candidatos — HFRS vs Control",
       subtitle = "Por tipo celular | * p<0.05  ** p<0.01  *** p<0.001 (Wilcoxon, BH)\nGris = no testado (n insuficiente)",
       x = "Tipo celular", y = "Receptor") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y   = element_text(size = 9),
        plot.subtitle = element_text(color = "grey40", size = 8),
        legend.position = "right")

all_plots[["20_Heatmap_receptores_DEG"]] <- p_heat_receptores

write.csv(
  deg_receptores_all %>% filter(DEG_class != "NS") %>%
    dplyr::select(gene, celltype, avg_log2FC, p_val_adj, DEG_class) %>%
    arrange(celltype, p_val_adj),
  paste0(results_dir, "Receptores_significativos_resumen.csv"), row.names = FALSE)

#================================================================
# BLOQUE 13: GUARDADO FINAL 
#================================================================
save_pub(paste0(output_dir, "Fig02_QC_violin.png"),
         all_plots[["02_QC_violin"]], width = 16, height = 6, dpi = DPI_STD)

save_pub(paste0(output_dir, "Fig03_QC_scatter.png"),
         all_plots[["03_QC_scatter"]], width = 14, height = 7, dpi = DPI_STD)

save_pub(paste0(output_dir, "Fig04b_PCA_pseudobulk.png"),
         all_plots[["04b_PCA_pseudobulk"]], width = 8, height = 7, dpi = DPI_HIGH)

save_pub(paste0(output_dir, "Fig05_UMAP_clusters.png"),
         all_plots[["05_UMAP_clusters"]], width = 10, height = 8, dpi = DPI_STD)

save_pub(paste0(output_dir, "Fig06_UMAP_SingleR.png"),
         all_plots[["06_UMAP_SingleR"]], width = 13, height = 9, dpi = DPI_STD)

save_pub(paste0(output_dir, "Fig07_Proporciones_clinicas.png"),
         all_plots[["07_Proporciones_clinicas"]], width = 8, height = 6, dpi = DPI_STD)

save_pub(paste0(output_dir, "Fig09_UMAP_inmune.png"),
         all_plots[["10_UMAP_inmune"]], width = 12, height = 9, dpi = DPI_STD)

if (!is.null(all_plots[["11_Feature_candidatos"]]))
  save_pub(paste0(output_dir, "Fig10_Feature_candidatos.png"),
           all_plots[["11_Feature_candidatos"]], width = 15, height = 14, dpi = DPI_STD)

if (!is.null(all_plots[["12_Dot_candidatos"]]))
  save_pub(paste0(output_dir, "Fig11_Dot_candidatos.png"),
           all_plots[["12_Dot_candidatos"]], width = 13, height = 8, dpi = DPI_STD)

if (!is.null(all_plots[["13_Feature_dirigidos"]]))
  save_pub(paste0(output_dir, "Fig12_Feature_dirigidos.png"),
           all_plots[["13_Feature_dirigidos"]], width = 14, height = 12, dpi = DPI_STD)

if (!is.null(all_plots[["14_Dot_dirigidos"]]))
  save_pub(paste0(output_dir, "Fig13_Dot_dirigidos.png"),
           all_plots[["14_Dot_dirigidos"]], width = 12, height = 8, dpi = DPI_STD)

if (!is.null(all_plots[["20_Heatmap_receptores_DEG"]])) {
  n_genes <- n_distinct(deg_receptores_all$gene)
  n_tipos <- n_distinct(deg_receptores_all$celltype)
  save_pub(paste0(output_dir, "Fig15_Heatmap_receptores_DEG.png"),
           all_plots[["20_Heatmap_receptores_DEG"]],
           width  = max(8, n_tipos * 1.5 + 3),
           height = max(6, n_genes * 0.5 + 2),
           dpi    = DPI_STD)
}
saveRDS(seurat_filtered, "data/processed/human/seurat_GSE161354_processed.rds")
write.csv(candidatos, paste0(results_dir, "candidatos_receptores_final.csv"), row.names = FALSE)

