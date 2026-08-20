# ==============================================================================
# PIPELINE bulk RNA-seq — GSE198751 (versión simplificada)
# Organismo: Myodes glareolus (topillo rojo — reservorio natural de PUUV)
# NOTA IMPORTANTE: los IDs de gen del conteo son ENSMUSG... (reads mapeados
# contra el genoma de ratón como proxy, ya que M. glareolus no tiene OrgDb
# propio). Toda la anotación/enriquecimiento de este script usa org.Mm.eg.db.
# ==============================================================================

rm(list = ls())

# ----------------------------
# 0. LIBRERÍAS Y CONFIGURACIÓN
# ----------------------------

library(DESeq2)
library(GEOquery)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)

raw_counts_path <- "data/raw/vector/GSE198751_RAW/GSE198751_raw_counts.txt"
output_dir      <- "figures/vector/GSE198751_Mglareolus/"
results_dir     <- "results/vector/GSE198751_Mglareolus/"

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(raw_counts_path), recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

DPI_STD  <- 300
DPI_HIGH <- 400

save_pub <- function(path, plot, width, height, dpi = DPI_STD) {
  path <- sub("\\.jpe?g$", ".png", path, ignore.case = TRUE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height,
                  units = "in", res = dpi, background = "white")
  } else {
    png(filename = path, width = width, height = height, units = "in",
        res = dpi, type = "cairo", bg = "white")
  }
  tryCatch(print(plot), error = function(e) {
    warning(sprintf("Error renderizando %s: %s", basename(path), e$message))
  })
  dev.off()
  cat(sprintf("  OK: %s\n", basename(path)))
}

# ------------------------------------------------------------------------------
# 1. DESCARGA Y LECTURA DE DATOS
# ------------------------------------------------------------------------------

if (!file.exists(raw_counts_path)) {
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE198751",
                "&format=file&file=GSE198751_raw_counts.txt.gz")
  download.file(url, destfile = raw_counts_path, mode = "wb")
}

data_raw <- read.table(raw_counts_path, header = TRUE, sep = "\t",
                       row.names = 1, check.names = FALSE)

data_counts <- data_raw %>%
  dplyr::select(-any_of(c("Chr", "Start", "End", "Strand", "Length")))

rownames(data_counts) <- sub("^gene:", "", rownames(data_counts))

cat(sprintf("Matriz: %d genes x %d muestras\n", nrow(data_counts), ncol(data_counts)))

coldata <- data.frame(
  condition = factor(sub("_[0-9]+$", "", colnames(data_counts)),
                     levels = c("NI", "PHV", "PUUV")),
  row.names = colnames(data_counts)
)
print(table(coldata$condition))

# ------------------------------------------------------------------------------
# 2. ANOTACIÓN DE GENES (org.Mm.eg.db — IDs son ENSMUSG)
# ------------------------------------------------------------------------------

annotation_file <- paste0(results_dir, "gene_annotation.rds")

if (!file.exists(annotation_file)) {
  gene_ids <- rownames(data_counts)
  annotation <- data.frame(
    gene_id = gene_ids,
    symbol  = mapIds(org.Mm.eg.db, keys = gene_ids, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first"),
    entrez  = mapIds(org.Mm.eg.db, keys = gene_ids, column = "ENTREZID",
                     keytype = "ENSEMBL", multiVals = "first"),
    stringsAsFactors = FALSE
  )
  saveRDS(annotation, annotation_file)
  cat(sprintf("Anotados: %d / %d genes\n",
              sum(!is.na(annotation$symbol) & annotation$symbol != ""),
              nrow(annotation)))
} 
annotation <- readRDS(annotation_file)

# ------------------------------------------------------------------------------
# 3. DESeq2, VST, PCA
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(data_counts)),
  colData   = coldata,
  design    = ~ condition
)

keep <- rowSums(counts(dds) >= 5) >= 3
dds  <- dds[keep, ]
cat(sprintf("Genes tras filtrado: %d\n", nrow(dds)))

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

pal_cond <- c("NI" = "#74C8EC", "PHV" = "#F5A623", "PUUV" = "#E8755A")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_data$sample <- rownames(pca_data)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.2, color = "black", max.overlaps = 20) +
  scale_color_manual(values = pal_cond, name = "Condición") +
  labs(title = "PCA — GSE198751 (Myodes glareolus)",
       subtitle = sprintf("PC1: %d%%  |  PC2: %d%%", percentVar[1], percentVar[2]),
       x = sprintf("PC1: %d%% varianza", percentVar[1]),
       y = sprintf("PC2: %d%% varianza", percentVar[2])) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

save_pub(paste0(output_dir, "Fig01_PCA.png"), p_pca, width = 8, height = 7)

# (heatmap de distancias eliminado — ver Fig04 para heatmap sobre DEGs)

# ------------------------------------------------------------------------------
# 4. CONTRASTES DESeq2: xtracción y anotación
# ------------------------------------------------------------------------------

contrastes <- list(
  PUUV_vs_NI  = c("condition", "PUUV", "NI"),
  PHV_vs_NI   = c("condition", "PHV",  "NI"),
  PHV_vs_PUUV = c("condition", "PHV",  "PUUV")
)

get_results <- function(contraste, nombre) {
  res_raw    <- results(dds, contrast = contraste, alpha = 0.05)
  res_shrink <- lfcShrink(dds, contrast = contraste, type = "ashr", res = res_raw)
  
  as.data.frame(res_shrink) %>%
    rownames_to_column("gene_id") %>%
    left_join(annotation, by = "gene_id") %>%
    filter(!is.na(padj)) %>%
    mutate(
      label_gene = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id),
      DEG_class  = factor(case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "NS"
      ), levels = c("Up", "Down", "NS")),
      contraste = nombre
    )
}

res_lista <- list()
for (nombre in names(contrastes)) {
  cat(sprintf("\nContraste: %s\n", nombre))
  res_df <- get_results(contrastes[[nombre]], nombre)
  cat(sprintf("  DEGs (padj<0.05, |LFC|>1): %d Up | %d Down\n",
              sum(res_df$DEG_class == "Up"), sum(res_df$DEG_class == "Down")))
  
  res_lista[[nombre]] <- res_df
  write.csv(res_df, paste0(results_dir, "DESeq2_", nombre, "_full.csv"), row.names = FALSE)
  write.csv(res_df %>% filter(DEG_class != "NS"),
            paste0(results_dir, "DESeq2_", nombre, "_DEGs.csv"), row.names = FALSE)
}

write.csv(bind_rows(res_lista), paste0(results_dir, "DESeq2_all_contrasts.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. VOLCANO PLOTS (genes etiquetados por significancia/magnitud)
# ------------------------------------------------------------------------------

contraste_meta <- list(
  PUUV_vs_NI  = list(titulo = "PUUV vs NI",  sub = "Virus patogénico en el reservorio",
                     col_up = "#d62728", col_down = "#1f77b4"),
  PHV_vs_NI   = list(titulo = "PHV vs NI",   sub = "Virus no patogénico en el reservorio",
                     col_up = "#E8755A", col_down = "#74C8EC"),
  PHV_vs_PUUV = list(titulo = "PHV vs PUUV", sub = "Patogénico vs no patogénico",
                     col_up = "#F5A623", col_down = "#7B5EA7")
)

# Función genérica de volcano. Si se pasa `highlight_genes`, esos genes se
# resaltan con un color propio y siempre se etiquetan (independientemente de
# si son DEG o no) — útil para paneles de interés (p. ej. mieloide/linfoide).
plot_volcano <- function(res_df, titulo, sub, col_up, col_down,
                         highlight_genes = NULL, highlight_name = "Interés",
                         highlight_color = "#2CA02C", n_label = 15) {
  
  df <- res_df %>% mutate(logp = -log10(padj + 1e-300))
  
  if (!is.null(highlight_genes)) {
    df <- df %>% mutate(
      grupo = case_when(
        symbol %in% highlight_genes ~ highlight_name,
        DEG_class == "Up"           ~ "Up",
        DEG_class == "Down"         ~ "Down",
        TRUE                        ~ "NS"
      ),
      grupo = factor(grupo, levels = c(highlight_name, "Up", "Down", "NS"))
    )
    colores  <- setNames(c(highlight_color, col_up, col_down, "grey75"),
                         c(highlight_name, "Up", "Down", "NS"))
    # Etiquetar: marcadores resaltados con señal + TODOS los DEGs significativos
    top_lab  <- bind_rows(
      df %>% filter(grupo == highlight_name) %>%
        filter(DEG_class != "NS" | abs(log2FoldChange) > 0.5),
      df %>% filter(DEG_class %in% c("Up", "Down"))   # DEGs siempre etiquetados
    ) %>% distinct(gene_id, .keep_all = TRUE)
  } else {
    df <- df %>% mutate(grupo = DEG_class)
    colores <- c("Up" = col_up, "Down" = col_down, "NS" = "grey75")
    # Etiquetar TODOS los DEGs significativos + top 10 por significancia global
    # No se limita a top 15 porque hay pocos DEGs
    top_lab <- bind_rows(
      df %>% filter(DEG_class != "NS"),          # todos los significativos
      df %>% slice_max(logp, n = 10)             # outliers extremos en Y
    ) %>% distinct(gene_id, .keep_all = TRUE)
  }
  
  ggplot(df, aes(x = log2FoldChange, y = logp, color = grupo)) +
    geom_point(data = df %>% filter(grupo %in% c("NS")), alpha = 0.25, size = 1.3) +
    geom_point(data = df %>% filter(!grupo %in% c("NS")), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = colores, drop = FALSE, name = NULL) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(
      data          = top_lab,
      aes(label     = label_gene),
      size          = 3.1,
      color         = "black",
      fontface      = "bold",
      max.overlaps  = Inf,        # nunca descartar etiquetas
      force         = 2,          
      force_pull    = 0.5,
      box.padding   = 0.35,
      segment.size  = 0.3,
      segment.color = "grey50",
      min.segment.length = 0.1
    ) +
    labs(title = titulo, subtitle = sub, x = "Log2 Fold Change",
         y = expression(-log[10](padj))) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", plot.subtitle = element_text(color = "grey40", size = 9))
}

plots_volcano <- list()
for (nombre in names(contrastes)) {
  m <- contraste_meta[[nombre]]
  p <- plot_volcano(res_lista[[nombre]], m$titulo, m$sub, m$col_up, m$col_down)
  plots_volcano[[nombre]] <- p
}

p_panel <- plots_volcano$PUUV_vs_NI /
  plots_volcano$PHV_vs_NI /
  plots_volcano$PHV_vs_PUUV +
  plot_annotation(
    title    = "Volcanos: GSE198751 (Myodes glareolus)",
    subtitle = "De arriba a abajo: PUUV vs NI, PHV vs NI, PHV vs PUUV\nUmbrales: |LFC| > 1, padj < 0.05 (DESeq2, ashr)",
    theme    = theme(
      plot.title    = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )
save_pub(paste0(output_dir, "Fig03_volcanos_panel.png"), p_panel, width = 12, height = 21, dpi = DPI_HIGH)

# ------------------------------------------------------------------------------
# 5b. VOLCANOS MIELOIDE / LINFOIDE — genes obtenidos automáticamente de MSigDB C8
#     (Mus musculus). Se usa el contraste principal PUUV vs NI.
#     MSigDB C8 contiene firmas de tipos celulares inmunes; se filtra por linaje.
# ------------------------------------------------------------------------------

immune_sets <- msigdbr(species = "Mus musculus", category = "C8") %>%
  filter(str_detect(gs_name, "MONOCYTE|NEUTROPHIL|MYELOID|BCELL|TCELL|NK_CELL|PLASMA"))

genes_mieloides_msig <- immune_sets %>%
  filter(str_detect(gs_name, "MONOCYTE|NEUTROPHIL|MYELOID")) %>%
  pull(gene_symbol) %>% unique()

genes_linfoides_msig <- immune_sets %>%
  filter(str_detect(gs_name, "BCELL|TCELL|NK_CELL|PLASMA")) %>%
  pull(gene_symbol) %>% unique()

cat(sprintf("MSigDB C8 — Mieloides: %d genes | Linfoides: %d genes\n",
            length(genes_mieloides_msig), length(genes_linfoides_msig)))

m0 <- contraste_meta[["PUUV_vs_NI"]]

p_mieloide <- plot_volcano(
  res_lista[["PUUV_vs_NI"]],
  titulo = "PUUV vs NI — linaje mieloide (MSigDB C8)",
  sub    = "Resaltados: genes de firmas mieloides (monocitos, neutrófilos, macrófagos)",
  col_up = m0$col_up, col_down = m0$col_down,
  highlight_genes = genes_mieloides_msig,
  highlight_name  = "Mieloide",
  highlight_color = "#8C4B2E"
)

p_linfoide <- plot_volcano(
  res_lista[["PUUV_vs_NI"]],
  titulo = "PUUV vs NI — linaje linfoide (MSigDB C8)",
  sub    = "Resaltados: genes de firmas linfoides (T, B, NK, plasma)",
  col_up = m0$col_up, col_down = m0$col_down,
  highlight_genes = genes_linfoides_msig,
  highlight_name  = "Linfoide",
  highlight_color = "#2E5C8C"
)

save_pub(paste0(output_dir, "Fig03b_volcano_mieloide.png"), p_mieloide, width = 10, height = 7, dpi = DPI_HIGH)
save_pub(paste0(output_dir, "Fig03c_volcano_linfoide.png"), p_linfoide, width = 10, height = 7, dpi = DPI_HIGH)

p_panel_ml <- p_mieloide | p_linfoide +
  plot_annotation(
    title    = "Genes mieloides y linfoides — PUUV vs NI (GSE198751)",
    subtitle = "Firmas de tipo celular: MSigDB C8 (Mus musculus)",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )
save_pub(paste0(output_dir, "Fig03d_panel_miel_linf.png"), p_panel_ml, width = 18, height = 7, dpi = DPI_HIGH)


# ------------------------------------------------------------------------------
# 6. HEATMAP DEGs — UNIÓN DE LOS 3 CONTRASTES
# Muestra todos los genes significativos en al menos un contraste.
# Las filas se anotan con en qué contraste(s) son DEG y en qué dirección.
# Las columnas se anotan por condición.
# ------------------------------------------------------------------------------

# Unión de DEGs de los 3 contrastes (significativos en al menos uno)
degs_union <- bind_rows(res_lista) %>%
  filter(DEG_class != "NS") %>%
  pull(gene_id) %>%
  unique()


if (length(degs_union) >= 5) {
  
  # Limitar a los top por significancia global si hay muchos
  N_MAX <- 80
  if (length(degs_union) > N_MAX) {
    degs_union <- bind_rows(res_lista) %>%
      filter(gene_id %in% degs_union) %>%
      group_by(gene_id) %>%
      summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
      arrange(min_padj) %>%
      slice_head(n = N_MAX) %>%
      pull(gene_id)
  }
  
  # Etiquetas de símbolo para filas
  gene_labels <- annotation %>%
    filter(gene_id %in% degs_union) %>%
    mutate(label = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id)) %>%
    dplyr::select(gene_id, label)
  
  mat_hm <- assay(vsd)[degs_union, ]
  rownames(mat_hm) <- make.unique(
    gene_labels$label[match(rownames(mat_hm), gene_labels$gene_id)]
  )
  mat_scaled <- t(scale(t(mat_hm)))
  
  # Anotación de columnas (condición)
  col_annot <- data.frame(
    Condition = coldata$condition,
    row.names = rownames(coldata)
  )
  
  # Anotación de filas: en qué contraste(s) es DEG y dirección
  row_annot <- lapply(names(contrastes), function(nombre) {
    res_lista[[nombre]] %>%
      filter(gene_id %in% degs_union) %>%
      dplyr::select(gene_id, DEG_class) %>%
      mutate(DEG_class = as.character(DEG_class)) %>%
      dplyr::rename(!!nombre := DEG_class)
  }) %>%
    Reduce(function(a, b) full_join(a, b, by = "gene_id"), .) %>%
    mutate(across(-gene_id, ~ replace_na(., "NS"))) %>%
    column_to_rownames("gene_id")
  
  # Alinear con filas de la matriz (por gene_id original antes de renombrar)
  row_annot_ordered <- row_annot[degs_union, , drop = FALSE]
  rownames(row_annot_ordered) <- rownames(mat_scaled)
  
  # Colores anotación filas
  deg_colors <- c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey90")
  ann_colors <- list(
    Condition   = pal_cond,
    PUUV_vs_NI  = deg_colors,
    PHV_vs_NI   = deg_colors,
    PHV_vs_PUUV = deg_colors
  )
  
  n_genes <- nrow(mat_scaled)
  pheatmap(
    mat_scaled,
    annotation_col    = col_annot,
    annotation_row    = row_annot_ordered,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = max(5, 9 - n_genes %/% 15),
    cluster_cols      = TRUE,
    cluster_rows      = TRUE,
    color = colorRampPalette(c("#313695", "white", "#a50026"))(100),
    main     = "DEGs — unión de los 3 contrastes (Myodes glareolus)",
    filename = paste0(output_dir, "Fig04_heatmap_DEGs_3contrastes.png"),
    width = 13, height = max(10, n_genes * 0.18 + 3)
  )
  
  
  # ------------------------------------------------------------------------------
  # 7. ENRIQUECIMIENTO FUNCIONAL (GO BP + KEGG, org.Mm.eg.db)
  # ------------------------------------------------------------------------------
  
  run_enrichment <- function(res_df, nombre) {
    
    genes_up   <- res_df %>% filter(DEG_class == "Up",   !is.na(entrez)) %>% pull(entrez) %>% unique()
    genes_down <- res_df %>% filter(DEG_class == "Down", !is.na(entrez)) %>% pull(entrez) %>% unique()
    universe   <- res_df %>% filter(!is.na(entrez)) %>% pull(entrez) %>% unique()
    
    cat(sprintf("  Universo: %d | UP: %d | DOWN: %d\n", length(universe), length(genes_up), length(genes_down)))
    
    out <- list()
    
    if (length(genes_up) >= 5) {
      out$ego_up <- tryCatch(
        enrichGO(genes_up, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 readable = TRUE) %>%
          clusterProfiler::simplify(cutoff = 0.4, by = "p.adjust", select_fun = min),
        error = function(e) NULL
      )
      if (!is.null(out$ego_up) && nrow(as.data.frame(out$ego_up)) > 0)
        write.csv(as.data.frame(out$ego_up), paste0(results_dir, "GO_BP_up_", nombre, ".csv"), row.names = FALSE)
    }
    
    if (length(genes_down) >= 5) {
      out$ego_down <- tryCatch(
        enrichGO(genes_down, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 readable = TRUE) %>%
          clusterProfiler::simplify(cutoff = 0.4, by = "p.adjust", select_fun = min),
        error = function(e) NULL
      )
      if (!is.null(out$ego_down) && nrow(as.data.frame(out$ego_down)) > 0)
        write.csv(as.data.frame(out$ego_down), paste0(results_dir, "GO_BP_down_", nombre, ".csv"), row.names = FALSE)
    }
    
    genes_all <- c(genes_up, genes_down)
    if (length(genes_all) >= 5) {
      out$ekegg <- tryCatch(
        enrichKEGG(genes_all, universe = universe, organism = "mmu",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05),
        error = function(e) NULL
      )
      if (!is.null(out$ekegg) && nrow(as.data.frame(out$ekegg)) > 0)
        write.csv(as.data.frame(out$ekegg), paste0(results_dir, "KEGG_", nombre, ".csv"), row.names = FALSE)
    }
    
    out
  }
  
  enrich_lista <- list()
  for (nombre in names(contrastes)) enrich_lista[[nombre]] <- run_enrichment(res_lista[[nombre]], nombre)
  
  # Dotplots de GO BP para cada contraste
  for (nombre in names(contrastes)) {
    e <- enrich_lista[[nombre]]
    if (!is.null(e$ego_up) && nrow(as.data.frame(e$ego_up)) > 0) {
      p <- dotplot(e$ego_up, showCategory = 15) +
        ggtitle(sprintf("GO BP UP — %s", nombre)) + theme_minimal(base_size = 10)
      save_pub(paste0(output_dir, "Fig05_GO_up_", nombre, ".png"), p, width = 10, height = 8)
    }
    if (!is.null(e$ego_down) && nrow(as.data.frame(e$ego_down)) > 0) {
      p <- dotplot(e$ego_down, showCategory = 15) +
        ggtitle(sprintf("GO BP DOWN — %s", nombre)) + theme_minimal(base_size = 10)
      save_pub(paste0(output_dir, "Fig05_GO_down_", nombre, ".png"), p, width = 10, height = 8)
    }
  }
  
  # ------------------------------------------------------------------------------
  # 7b. GSEA — PUUV vs NI (complemento al ORA cuando hay pocos DEGs)
  #     Con solo ~28 DEGs significativos, el ORA tiene poco poder.
  #     GSEA usa el ranking completo de todos los genes → más robusto con n bajo.
  # ------------------------------------------------------------------------------
  
  gene_rank_PUUV <- res_lista[["PUUV_vs_NI"]] %>%
    filter(!is.na(entrez), !is.na(log2FoldChange), !is.na(pvalue)) %>%
    mutate(rank_score = sign(log2FoldChange) * (-log10(pvalue + 1e-300))) %>%
    arrange(desc(rank_score)) %>%
    dplyr::select(entrez, rank_score) %>%
    deframe()
  gene_rank_PUUV <- gene_rank_PUUV[!duplicated(names(gene_rank_PUUV))]
  
  if (length(gene_rank_PUUV) >= 100) {
    gsea_PUUV <- tryCatch(
      gseGO(geneList      = gene_rank_PUUV,
            OrgDb         = org.Mm.eg.db,
            keyType       = "ENTREZID",
            ont           = "BP",
            minGSSize     = 15,
            maxGSSize     = 500,
            pvalueCutoff  = 0.05,
            pAdjustMethod = "BH",
            verbose       = FALSE),
      error = function(e) { cat(sprintf("  GSEA falló: %s\n", e$message)); NULL }
    )
    
    if (!is.null(gsea_PUUV) && nrow(as.data.frame(gsea_PUUV)) > 0) {
      write.csv(as.data.frame(gsea_PUUV),
                paste0(results_dir, "GSEA_GO_BP_PUUV_vs_NI.csv"), row.names = FALSE)
      p_gsea <- dotplot(gsea_PUUV, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle("GSEA GO BP — PUUV vs NI (Myodes glareolus)") +
        theme_minimal(base_size = 10) +
        theme(axis.text.y = element_text(size = 8))
      save_pub(paste0(output_dir, "Fig05b_GSEA_PUUV_vs_NI.png"), p_gsea, width = 14, height = 10)
    }
  }
  
  # -----------------
  # 8. GUARDADO FINAL
  # -----------------
  write.csv(as.data.frame(counts(dds, normalized = TRUE)),
            paste0(results_dir, "normalized_counts.csv"))
  
  resumen_degs <- lapply(names(contrastes), function(nombre) {
    res_lista[[nombre]] %>% count(DEG_class) %>% mutate(contraste = nombre)
  }) %>% bind_rows() %>% pivot_wider(names_from = DEG_class, values_from = n, values_fill = 0)
  
  write.csv(resumen_degs, paste0(results_dir, "Resumen_DEGs_por_contraste.csv"), row.names = FALSE)
  saveRDS(dds, "data/processed/dds_GSE198751.rds")