# ==============================================================================
# PIPELINE bulk RNA-seq — PRJNA258076 
# Organismo: Oligoryzomys longicaudatus 
# Tejido: Bazo, RNA total, Ribo-Zero, HiSeq 2000, paired-end 2x100nt
#
# Diseño experimental:
#   ANDV (infectados):
#     RR18 —  carga viral baja -> PERSISTENTE
#     RR29 —carga viral alta   -> AGUDO
#     RR31 —  carga viral baja -> PERSISTENTE
# Controles
#     RR20 
#     RR30
#
# Contrastes (replicando el artículo original):
#   1. Persistente: RR18+RR31 vs RR20+RR30  (n=2 vs n=2)
#   2. Agudo:       RR29       vs RR20+RR30  (n=1 vs n=2)
#
# Datos de entrada:
#   counts estimados por Salmon sobre Mus musculus Ensembl
#   (alternativa rápida a Trinity de novo del artículo original)
#
# PARÁMETROS ORIGINALES replicados en DESeq2:
#   Persistente: method="pooled", sharingMode="maximum", fitType="local"
#   Agudo (RR29): method="blind", sharingMode="fit-only", fitType="local"
#
# Filtro mínimo de counts: sum >= 200 reads (igual que el original)
# Umbral DEG: FDR (padj) < 0.05 (igual que el original)
# ==============================================================================


rm(list = ls())
# ------------------------------------------------------------------------------
# 0. LIBRERÍAS Y CONFIGURACIÓN
# ------------------------------------------------------------------------------

library(DESeq2)
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
library(RColorBrewer)
library(GenomicFeatures)
library(tximport)

raw_counts_dir <- "data/raw/vector/PRJNA258076_RAW/"
output_dir     <- "figures/vector/PRJNA258076/"
results_dir    <- "results/vector/PRJNA258076/"

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

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

write_csv_safe <- function(df, path, ...) {
  culpables <- names(df)[sapply(df, is.list)]
  if (length(culpables) > 0) {
    cat(sprintf("  Aviso: columnas tipo lista detectadas en %s -> %s (aplanando)\n",
                basename(path), paste(culpables, collapse = ", ")))
    df <- df %>%
      mutate(across(where(is.list), ~ vapply(
        ., function(x) if (length(x) == 0 || all(is.na(x))) NA_character_ else paste(x, collapse = ";"),
        character(1)
      )))
  }
  write.csv(df, path, ...)
}

# ------------------------------------------------------------------------------
# 1. CARGA DE CUANTIFICACIONES SALMON (tximport)
# ------------------------------------------------------------------------------

muestras <- c("RR18","RR20","RR29","RR30","RR31")

files <- file.path(
  "data/processed/salmon_quant",
  muestras,
  "quant.sf"
)

names(files) <- muestras

if (any(!file.exists(files))) {
  stop("Faltan archivos quant.sf")
}

gtf <- "data/reference/Mus_musculus.GRCm39.113.gtf"

txdb <- makeTxDbFromGFF(gtf)

k <- keys(txdb, keytype="TXNAME")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys=k,
  keytype="TXNAME",
  columns="GENEID"
)

colnames(tx2gene) <- c("TXNAME","GENEID")

txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene,
  ignoreTxVersion=TRUE
)

count_matrix <- round(txi$counts)

cat(sprintf(
  "Matriz: %d genes x %d muestras\n",
  nrow(count_matrix),
  ncol(count_matrix)
))

filtro_original <- rowSums(count_matrix) >= 200
count_matrix <- count_matrix[filtro_original,]

cat(sprintf(
  "Tras filtro >=200 reads: %d genes\n",
  nrow(count_matrix)
))

# ------------------------------------------------------------------------------
# 2. ANOTACIÓN
# ------------------------------------------------------------------------------

gene_ids <- rownames(count_matrix)

annotation <- data.frame(
  gene_id = gene_ids,
  stringsAsFactors = FALSE
)

annotation$symbol <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys=gene_ids,
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first"
)

annotation$entrez <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys=gene_ids,
  column="ENTREZID",
  keytype="ENSEMBL",
  multiVals="first"
)

annotation$symbol[is.na(annotation$symbol)] <- annotation$gene_id

# ------------------------------------------------------------------------------
# 3. METADATOS
# ------------------------------------------------------------------------------

coldata <- data.frame(
  muestra   = muestras,
  estado    = factor(c("Persistente", "Control", "Agudo", "Control", "Persistente"),
                     levels = c("Control", "Persistente", "Agudo")),
  sexo      = c("F", "M", "M", "F", "M"),
  row.names = muestras,
  stringsAsFactors = FALSE
)

pal_estado <- c("Control" = "#74C8EC", "Persistente" = "#E8755A", "Agudo" = "#B0304A")

# ------------------------------------------------------------------------------
# 4. PCA
# ------------------------------------------------------------------------------

dds_qc <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ estado
)
dds_qc <- estimateSizeFactors(dds_qc)
vsd_qc <- varianceStabilizingTransformation(dds_qc, blind = TRUE)

pca_data   <- plotPCA(vsd_qc, intgroup = "estado", returnData = TRUE)
pca_data$muestra <- rownames(pca_data)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = estado, label = muestra)) +
  geom_point(size = 7, alpha = 0.9) +
  geom_text_repel(size = 3.8, color = "black", fontface = "bold",
                  box.padding = 0.5, max.overlaps = 20) +
  scale_color_manual(values = pal_estado, name = "Estado clínico") +
  annotate("text", x = min(pca_data$PC1), y = max(pca_data$PC2),
           label = "n=5 individuos silvestres outbred",
           hjust = 0, vjust = 1, size = 3, color = "grey40", fontface = "italic") +
  labs(
    title    = "PCA — PRJNA258076 (O. longicaudatus, bazo)",
    subtitle = sprintf("PC1: %d%%  |  PC2: %d%%", percentVar[1], percentVar[2]),
    x = sprintf("PC1: %d%%", percentVar[1]),
    y = sprintf("PC2: %d%%", percentVar[2])
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", plot.subtitle = element_text(color = "grey40", size = 9))

save_pub(paste0(output_dir, "Fig01_PCA.png"), p_pca, width = 9, height = 7, dpi = DPI_HIGH)

# ------------------------------------------------------------------------------
# 5. ANÁLISIS DIFERENCIAL
# ------------------------------------------------------------------------------

# CONTRASTE 1: PERSISTENTE
cat("\n=== CONTRASTE PERSISTENTE ===\n")
muestras_pers <- c("RR18", "RR20", "RR30", "RR31")
counts_pers   <- count_matrix[, muestras_pers]
coldata_pers  <- coldata[muestras_pers, , drop = FALSE]
coldata_pers$grupo <- factor(
  ifelse(coldata_pers$muestra %in% c("RR18", "RR31"), "Persistente", "Control"),
  levels = c("Control", "Persistente")
)

dds_pers <- DESeqDataSetFromMatrix(
  countData = counts_pers, colData = coldata_pers, design = ~ grupo
)
dds_pers <- estimateSizeFactors(dds_pers)
dds_pers <- estimateDispersions(dds_pers, fitType = "local")
dds_pers <- nbinomWaldTest(dds_pers)

res_pers_raw <- results(dds_pers,
                        contrast  = c("grupo", "Persistente", "Control"),
                        alpha     = 0.05)

res_pers_df <- as.data.frame(res_pers_raw) %>%
  rownames_to_column("gene_id") %>%
  left_join(annotation, by = "gene_id") %>%
  filter(!is.na(padj)) %>%
  mutate(
    across(everything(), ~ if (is.list(.)) sapply(., paste, collapse = ";") else .),
    label_gene = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id),
    DEG_class  = factor(case_when(
      padj < 0.05 & log2FoldChange > 0  ~ "Up",
      padj < 0.05 & log2FoldChange < 0  ~ "Down",
      TRUE ~ "NS"
    ), levels = c("Up", "Down", "NS")),
    contraste = "Persistente"
  )

n_up_pers   <- sum(res_pers_df$DEG_class == "Up",   na.rm = TRUE)
n_down_pers <- sum(res_pers_df$DEG_class == "Down",  na.rm = TRUE)
cat(sprintf("DEGs (padj<0.05): %d Up | %d Down\n", n_up_pers, n_down_pers))

write_csv_safe(res_pers_df,
               paste0(results_dir, "DESeq2_Persistente_vs_Control_full.csv"), row.names = FALSE)
write_csv_safe(res_pers_df %>% filter(DEG_class != "NS"),
               paste0(results_dir, "DESeq2_Persistente_vs_Control_DEGs.csv"), row.names = FALSE)

# CONTRASTE 2: AGUDO
cat("\n=== CONTRASTE AGUDO ===\n")
muestras_agudo <- c("RR29", "RR20", "RR30")
counts_agudo   <- count_matrix[, muestras_agudo]
coldata_agudo  <- coldata[muestras_agudo, , drop = FALSE]
coldata_agudo$grupo <- factor(
  ifelse(coldata_agudo$muestra == "RR29", "Agudo", "Control"),
  levels = c("Control", "Agudo")
)

dds_agudo <- DESeqDataSetFromMatrix(
  countData = counts_agudo, colData = coldata_agudo, design = ~ grupo
)
dds_agudo <- estimateSizeFactors(dds_agudo)
dds_agudo <- estimateDispersions(dds_agudo)
dds_agudo <- nbinomWaldTest(dds_agudo)

res_agudo_raw <- results(dds_agudo,
                         contrast = c("grupo", "Agudo", "Control"),
                         alpha    = 0.05)

res_agudo_df <- as.data.frame(res_agudo_raw) %>%
  rownames_to_column("gene_id") %>%
  left_join(annotation, by = "gene_id") %>%
  filter(!is.na(padj)) %>%
  mutate(
    across(everything(), ~ if (is.list(.)) sapply(., paste, collapse = ";") else .),
    label_gene = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id),
    DEG_class  = factor(case_when(
      padj < 0.05 & log2FoldChange > 0  ~ "Up",
      padj < 0.05 & log2FoldChange < 0  ~ "Down",
      TRUE ~ "NS"
    ), levels = c("Up", "Down", "NS")),
    contraste = "Agudo"
  )

n_up_agudo   <- sum(res_agudo_df$DEG_class == "Up",   na.rm = TRUE)
n_down_agudo <- sum(res_agudo_df$DEG_class == "Down",  na.rm = TRUE)
cat(sprintf("DEGs (padj<0.05): %d Up | %d Down\n", n_up_agudo, n_down_agudo))

write_csv_safe(res_agudo_df,
               paste0(results_dir, "DESeq2_Agudo_vs_Control_full.csv"), row.names = FALSE)
write_csv_safe(res_agudo_df %>% filter(DEG_class != "NS"),
               paste0(results_dir, "DESeq2_Agudo_vs_Control_DEGs.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. VOLCANOS
# ------------------------------------------------------------------------------

genes_articulo <- c("Casp1", "Stat2", "Ifih1", "Ddx58", "Ifit1", "Ifit2", "Mx1", 
                    "Mx2", "Isg15", "Oas1", "Rsad2", "Irf7", "Irf3", "Sting1",
                    "Il1b", "Tnf", "Ifng", "Cd55", "Tlr3", "Tlr7", "Tlr8")

plot_volcano_prjna <- function(res_df, titulo, sub,
                               col_up = "#E8755A", col_down = "#4C72B0") {
  df <- res_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    mutate(logp = -log10(padj + 1e-300))
  
  lab_up   <- df %>% filter(DEG_class == "Up") %>% arrange(padj) %>% slice_head(n = 15)
  lab_down <- df %>% filter(DEG_class == "Down") %>% arrange(padj) %>% slice_head(n = 10)
  lab_art  <- df %>% filter(symbol %in% genes_articulo, DEG_class != "NS")
  top_lab  <- bind_rows(lab_up, lab_down, lab_art) %>% distinct(gene_id, .keep_all = TRUE)
  
  df <- df %>% mutate(
    label_plot  = ifelse(gene_id %in% top_lab$gene_id, label_gene, ""),
    es_articulo = symbol %in% genes_articulo & DEG_class != "NS"
  )
  
  ggplot(df, aes(x = log2FoldChange, y = logp, color = DEG_class)) +
    geom_point(data = df %>% filter(DEG_class == "NS"),
               alpha = 0.2, size = 0.9, color = "grey75") +
    geom_point(data = df %>% filter(DEG_class != "NS" & !es_articulo),
               alpha = 0.75, size = 1.6) +
    geom_point(data = df %>% filter(es_articulo),
               size = 3, shape = 21, stroke = 1.2,
               aes(fill = DEG_class), color = "black") +
    scale_color_manual(values = c("Up" = col_up, "Down" = col_down, "NS" = "grey75"),
                       drop = FALSE, name = NULL) +
    scale_fill_manual(values  = c("Up" = col_up, "Down" = col_down, "NS" = "grey75"),
                      drop = FALSE, guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(
      aes(label = label_plot), size = 2.8,
      color = "black", fontface = "bold",
      max.overlaps = 30, force = 3, force_pull = 0.3,
      box.padding = 0.4, segment.size = 0.25,
      segment.color = "grey50", min.segment.length = 0.1, seed = 42
    ) +
    labs(title = titulo, subtitle = sub,
         x = "Log2 Fold Change", y = expression(-log[10](padj))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.subtitle = element_text(color = "grey40", size = 8))
}

p_volc_pers <- plot_volcano_prjna(
  res_pers_df,
  titulo = "Persistente (RR18+RR31) vs Control",
  sub    = sprintf("n=%d Up | n=%d Down", n_up_pers, n_down_pers),
  col_up = "#E8755A", col_down = "#4C72B0"
)

p_volc_agudo <- plot_volcano_prjna(
  res_agudo_df,
  titulo = "Agudo (RR29) vs Control",
  sub    = sprintf("n=%d Up | n=%d Down", n_up_agudo, n_down_agudo),
  col_up = "#B0304A", col_down = "#4C72B0"
)

p_panel <- (p_volc_pers | p_volc_agudo) +
  plot_annotation(
    title = "Expresión diferencial — O. longicaudatus + ANDV (PRJNA258076)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
save_pub(paste0(output_dir, "Fig04_Volcano_Panel.png"), p_panel, width = 18, height = 7)

# ------------------------------------------------------------------------------
# 7. ENRIQUECIMIENTO FUNCIONAL
# ------------------------------------------------------------------------------

run_enrichment_prjna <- function(res_df, nombre) {
  cat(sprintf("\nEnriquecimiento: %s\n", nombre))
  
  universe <- res_df %>% filter(!is.na(entrez)) %>% pull(entrez) %>% unique()
  genes_up <- res_df %>% filter(DEG_class == "Up", !is.na(entrez)) %>% pull(entrez) %>% unique()
  genes_dn <- res_df %>% filter(DEG_class == "Down", !is.na(entrez)) %>% pull(entrez) %>% unique()
  
  cat(sprintf("  Universo: %d | UP: %d | DOWN: %d\n",
              length(universe), length(genes_up), length(genes_dn)))
  
  out <- list()
  
  if (length(genes_up) >= 3) {
    out$ego_up <- tryCatch(
      enrichGO(genes_up, universe = universe, OrgDb = org.Mm.eg.db,
               keyType = "ENTREZID", ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.1,
               qvalueCutoff  = 0.3, readable = TRUE) %>%
        simplify(cutoff = 0.5, by = "p.adjust", select_fun = min),
      error = function(e) NULL
    )
  }
  
  if (length(genes_dn) >= 3) {
    out$ego_dn <- tryCatch(
      enrichGO(genes_dn, universe = universe, OrgDb = org.Mm.eg.db,
               keyType = "ENTREZID", ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.1,
               qvalueCutoff  = 0.3, readable = TRUE) %>%
        simplify(cutoff = 0.5, by = "p.adjust", select_fun = min),
      error = function(e) NULL
    )
  }
  
  for (tipo in c("ego_up", "ego_dn")) {
    if (!is.null(out[[tipo]]) && nrow(as.data.frame(out[[tipo]])) > 0) {
      write.csv(as.data.frame(out[[tipo]]),
                paste0(results_dir, tipo, "_", nombre, ".csv"), row.names = FALSE)
    }
  }
  
  out
}

enrich_pers  <- run_enrichment_prjna(res_pers_df,  "Persistente")
enrich_agudo <- run_enrichment_prjna(res_agudo_df, "Agudo")

# ------------------------------------------------------------------------------
# 7b. FIGURA ÚNICA COMPARATIVA — GO BP 
#     (Up/Down combinados en un mismo eje, Agudo vs Persistente lado a lado)
# ------------------------------------------------------------------------------

# Clasificación funcional simple por palabras clave en la descripción GO.
# Ajustar/ampliar esta lista si aparecen términos que caen todos en "Otros".
clasificar_categoria <- function(term) {
  term_low <- tolower(term)
  case_when(
    str_detect(term_low, "interferon|antiviral|defense response to virus|response to virus|viral") ~ "Respuesta antiviral / IFN",
    str_detect(term_low, "inflamm|cytokine|chemokine|acute phase|tumor necrosis|interleukin") ~ "Inflamación / citoquinas",
    str_detect(term_low, "t cell|b cell|lymphocyte|immune response|antigen|leukocyte|natural killer") ~ "Regulación inmune",
    str_detect(term_low, "coagulat|hemostasis|platelet|wound heal|vascular|angiogenesis|endothel") ~ "Vascular / hemostasia",
    str_detect(term_low, "cell cycle|mitotic|dna replication|chromosome|chromatin|centromere|meiotic") ~ "Ciclo celular / DNA",
    str_detect(term_low, "migration|chemotaxis|adhesion|locomotion|mesenchyme") ~ "Migración / adhesión",
    str_detect(term_low, "apoptosis|apoptotic|programmed cell death") ~ "Apoptosis",
    str_detect(term_low, "translation|ribosom") ~ "Traducción / ribosoma",
    str_detect(term_low, "mapk|signaling pathway|signal transduction") ~ "Señalización",
    TRUE ~ "Otros"
  )
}

pal_categoria <- c(
  "Respuesta antiviral / IFN" = "#2E86AB",
  "Inflamación / citoquinas"  = "#D4A017",
  "Regulación inmune"         = "#7B68C8",
  "Vascular / hemostasia"     = "#B0304A",
  "Ciclo celular / DNA"       = "#B5651D",
  "Migración / adhesión"      = "#1E8A6E",
  "Apoptosis"                 = "#8B4789",
  "Traducción / ribosoma"     = "#4C72B0",
  "Señalización"              = "#5F9E7C",
  "Otros"                     = "grey50"
)

N_TERMINOS_LOLLIPOP <- 10  # top N términos por dirección (up/down) y por contraste

extraer_df_go <- function(ego_obj, direccion, nombre, n_top = N_TERMINOS_LOLLIPOP) {
  if (is.null(ego_obj)) return(NULL)
  df <- as.data.frame(ego_obj)
  if (nrow(df) == 0) return(NULL)
  
  df %>%
    mutate(
      GeneRatio_num    = sapply(GeneRatio, function(x) eval(parse(text = x))),
      GeneRatio_signed = if (direccion == "Down") -GeneRatio_num else GeneRatio_num,
      categoria        = clasificar_categoria(Description),
      contraste        = nombre,
      direccion        = direccion
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top)
}

df_go_todo <- bind_rows(
  extraer_df_go(enrich_agudo$ego_up, "Up",   "Agudo"),
  extraer_df_go(enrich_agudo$ego_dn, "Down", "Agudo"),
  extraer_df_go(enrich_pers$ego_up,  "Up",   "Persistente"),
  extraer_df_go(enrich_pers$ego_dn,  "Down", "Persistente")
)

if (!is.null(df_go_todo) && nrow(df_go_todo) > 0) {
  # Orden compartido de términos en el eje Y: por categoría y luego por GeneRatio medio
  orden_terminos <- df_go_todo %>%
    group_by(Description, categoria) %>%
    summarise(orden_val = mean(GeneRatio_signed), .groups = "drop") %>%
    arrange(categoria, orden_val) %>%
    pull(Description)
  
  df_go_todo <- df_go_todo %>%
    mutate(
      Description = factor(Description, levels = orden_terminos),
      contraste   = factor(contraste, levels = c("Agudo", "Persistente"),
                           labels = c("Agudo\nvs Control", "Persistente\nvs Control"))
    )
  
  lim_x <- max(abs(df_go_todo$GeneRatio_signed)) * 1.15
  
  p_go_comparativo <- ggplot(df_go_todo, aes(x = GeneRatio_signed, y = Description)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey20") +
    geom_segment(aes(x = 0, xend = GeneRatio_signed, y = Description, yend = Description,
                     color = categoria), linewidth = 0.7) +
    geom_point(aes(size = Count, color = categoria), alpha = 0.9) +
    scale_color_manual(values = pal_categoria, name = "Categoría funcional") +
    scale_size_continuous(name = "Count", range = c(2, 8)) +
    scale_x_continuous(limits = c(-lim_x, lim_x), labels = abs) +
    facet_wrap(~ contraste, nrow = 1) +
    labs(
      title = "GO Biological Process — Agudo vs Control  |  Persistente vs Control",
      subtitle = sprintf("PRJNA258076 (O. longicaudatus + ANDV) | top %d términos por dirección | padj<0.1 | simplify cutoff=0.5",
                         N_TERMINOS_LOLLIPOP),
      x = "GeneRatio     \u2190  Down            Up  \u2192",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text     = element_text(face = "bold", size = 12),
      panel.spacing  = unit(1.5, "lines"),
      plot.subtitle  = element_text(color = "grey40", size = 9),
      legend.position = "right"
    )
  
  altura_fig <- max(6, 0.35 * length(orden_terminos) + 2.5)
  
  save_pub(paste0(output_dir, "Fig05_GO_comparativo.png"), p_go_comparativo,
           width = 14, height = altura_fig)
} else {
  cat("  Sin términos GO disponibles para generar la figura comparativa.\n")
}

# RESUMEN FINAL

cat(sprintf("Persistente: %d Up | %d Down\n", n_up_pers, n_down_pers))
cat(sprintf("Agudo: %d Up | %d Down\n", n_up_agudo, n_down_agudo))
cat(sprintf("\nFiguras: %s\n", output_dir))
cat(sprintf("Resultados: %s\n", results_dir))