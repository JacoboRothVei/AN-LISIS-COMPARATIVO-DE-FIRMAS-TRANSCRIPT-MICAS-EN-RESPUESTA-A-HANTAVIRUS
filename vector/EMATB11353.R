# ==============================================================================
# PIPELINE bulk RNA-seq — E-MTAB-11353
# Organismo: Mus musculus | mBMDM (macrófagos primarios de médula ósea)
# Estudio: RNA-seq de mBMDM infectados con HTNV (Hantaan virus, MOI=1) vs
#          mock (virus inactivado con Co60), time-course 0/12/24/36 hpi
#   A (0h)  = Control: HTNV inactivado con Co60 (mock, no infeccioso)
#   B (12h) = HTNV vivo, 12h post-infección
#   C (24h) = HTNV vivo, 24h post-infección
#   D (36h) = HTNV vivo, 36h post-infección
# Contrastes: cada tiempo de infección vs control [A(0h)]
#
# Los archivos son a nivel de TRANSCRITO (una fila por
# transcrito, columna "Gene" con el gen al que pertenece, columna "{sample}_COUNT"
# con el count). La columna "Gene" mezcla IDs Ensembl (ENSMUSG...) y símbolos ya
# resueltos (RefSeq/otros), el pipeline agrega a nivel de gen y anota ambos casos.
# ==============================================================================

rm(list = ls())
file.remove("results/vector/EMTAB11353/gene_annotation.rds")

# -----------------------------
# 0. LIBRERÍAS Y CONFIGURACIÓN
# -----------------------------

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
library(msigdbr)
library(rrvgo)    # rrvgo agrupa términos GO redundantes por similitud semántica

# Paths
raw_counts_dir <- "data/raw/vector/EMTAB11353_RAW/"
output_dir     <- "figures/vector/EMTAB11353_Mmusculus/"
results_dir    <- "results/vector/EMTAB11353_Mmusculus/"

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

#Caldiad de las figuras
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
  print(plot)
  dev.off()
}

# -------------------------------------------------------------------
# 1. LECTURA Y AGREGACIÓN DE COUNTS (transcrito -> gen, por muestra)
# -------------------------------------------------------------------
sample_files <- list.files(raw_counts_dir, pattern = "^Transcript_expression-.*\\.txt$",
                           full.names = TRUE)
stopifnot(length(sample_files) > 0)

sample_names <- sub("^Transcript_expression-(.*)\\.txt$", "\\1", basename(sample_files))
cat(sprintf("Muestras detectadas (%d): %s\n", length(sample_names), paste(sample_names, collapse = ", ")))

# Suma los counts de todos los transcritos de un mismo gen (columna "Gene").
# La columna "Gene" mezcla IDs Ensembl con versión (ENSMUSGxxxxxxxxxxx.N) y
# símbolos ya resueltos (RefSeq/otros) 
read_gene_counts <- function(filepath) {
  df <- read.table(filepath, header = TRUE, sep = "\t", quote = "",
                   comment.char = "", stringsAsFactors = FALSE)
  count_col <- grep("_COUNT$", colnames(df), value = TRUE)[1]
  
  gene_clean <- ifelse(
    grepl("^ENSMUSG[0-9]+\\.[0-9]+$", df$Gene),
    sub("\\.[0-9]+$", "", df$Gene),
    df$Gene
  )
  
  agg <- aggregate(df[[count_col]], by = list(gene_id = gene_clean), FUN = sum)
  colnames(agg) <- c("gene_id", "count")
  agg
}

#Agregando counts a nivel de gen por muestra
lista_counts <- lapply(sample_files, read_gene_counts)
names(lista_counts) <- sample_names

data_counts <- Reduce(function(a, b) full_join(a, b, by = "gene_id"),
                      Map(function(df, nm) setNames(df, c("gene_id", nm)),
                          lista_counts, names(lista_counts))) %>%
  mutate(across(-gene_id, ~ replace_na(., 0))) %>%
  column_to_rownames("gene_id")

data_counts <- data_counts[, sample_names]  # asegurar orden

cat(sprintf("Matriz final: %d genes x %d muestras\n", nrow(data_counts), ncol(data_counts)))

# ------------------------------------------------------------------------------
# 2. METADATOS (coldata) — a partir del prefijo de letra de cada muestra
#    (A=0h control mock, B=12h, C=24h, D=36h HTNV)
# ------------------------------------------------------------------------------

condition_map <- c(A = "Control_0h", B = "HTNV_12h", C = "HTNV_24h", D = "HTNV_36h")

coldata <- data.frame(
  condition = factor(
    condition_map[substr(colnames(data_counts), 1, 1)],
    levels = c("Control_0h", "HTNV_12h", "HTNV_24h", "HTNV_36h")
  ),
  row.names = colnames(data_counts)
)
print(table(coldata$condition))
stopifnot(!any(is.na(coldata$condition)))

# ------------------------------------------------------------------------------
# 3. ANOTACIÓN DE GENES (org.Mm.eg.db maneja IDs mixtos Ensembl/símbolo)
# ------------------------------------------------------------------------------

annotation_file <- paste0(results_dir, "gene_annotation.rds")

if (!file.exists(annotation_file)) {
  gene_ids  <- rownames(data_counts)
  is_ensembl <- grepl("^ENSMUSG[0-9]+$", gene_ids)
  
  cat(sprintf("IDs Ensembl: %d | símbolos/otros: %d\n", sum(is_ensembl), sum(!is_ensembl)))
  
  symbol_resuelto <- gene_ids
  symbol_resuelto[is_ensembl] <- mapIds(org.Mm.eg.db, keys = gene_ids[is_ensembl],
                                        column = "SYMBOL", keytype = "ENSEMBL",
                                        multiVals = "first")
  
  annotation <- data.frame(
    gene_id = gene_ids,
    symbol  = symbol_resuelto,
    stringsAsFactors = FALSE
  )
  
  # Entrez a partir del símbolo ya resuelto (uniforme para ambos orígenes de ID)
  annotation$entrez <- mapIds(org.Mm.eg.db, keys = annotation$symbol,
                              column = "ENTREZID", keytype = "SYMBOL",
                              multiVals = "first")
  
  saveRDS(annotation, annotation_file)
  cat(sprintf("Anotados: %d / %d genes (símbolo no NA)\n",
              sum(!is.na(annotation$symbol) & annotation$symbol != ""),
              nrow(annotation)))
} else {
  cat("Cargando anotación desde caché.\n")
}
annotation <- readRDS(annotation_file)


if (any(sapply(annotation, is.list))) {
  cat("  Aviso: anotación en caché tenía columnas tipo lista; aplanando...\n")
  annotation <- annotation %>%
    mutate(across(where(is.list), ~ vapply(
      ., function(x) if (length(x) == 0 || all(is.na(x))) NA_character_ else paste(x, collapse = ";"),
      character(1)
    )))
}

# ------------------------------------------------------------------------------
# 4. DESeq2, VST (estabuliza relación media-varianza), PCA
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

pal_cond <- c("Control_0h" = "#74C8EC", "HTNV_12h" = "#F5A623",
              "HTNV_24h" = "#E8755A", "HTNV_36h" = "#B0304A")

#PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_data$sample <- rownames(pca_data)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.2, color = "black", max.overlaps = 20) +
  scale_color_manual(values = pal_cond, name = "Condición") +
  labs(title = "PCA — E-MTAB-11353 (mBMDM, time-course HTNV)",
       subtitle = sprintf("PC1: %d%%  |  PC2: %d%%", percentVar[1], percentVar[2]),
       x = sprintf("PC1: %d%% varianza", percentVar[1]),
       y = sprintf("PC2: %d%% varianza", percentVar[2])) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

save_pub(paste0(output_dir, "Fig01_PCA.png"), p_pca, width = 8, height = 7)

# ------------------------------------------------------------------------------
# 5. CONTRASTES DESeq2 — cada tiempo de infección vs control
# ------------------------------------------------------------------------------

contrastes <- list(
  HTNV_12h_vs_Control = c("condition", "HTNV_12h", "Control_0h"),
  HTNV_24h_vs_Control = c("condition", "HTNV_24h", "Control_0h"),
  HTNV_36h_vs_Control = c("condition", "HTNV_36h", "Control_0h")
)


# Contraste "principal" de referencia (el más tardío/completo), usado en el
# volcano global, GSEA y paneles mieloide/linfoide. 
nombre_principal <- "HTNV_36h_vs_Control"

get_results <- function(contraste, nombre) {
  res_raw    <- results(dds, contrast = contraste, alpha = 0.05)
  res_shrink <- lfcShrink(dds, contrast = contraste, type = "ashr", res = res_raw)
  
  out <- as.data.frame(res_shrink) %>%
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
  
  # aplanar cualquier columna tipo lista que pudiera colarse, para que write.csv()
  # nunca falle más adelante con "tipo no implementado 'list' en 'EncodeElement'".
  if (any(sapply(out, is.list))) {
    out <- out %>%
      mutate(across(where(is.list), ~ vapply(
        ., function(x) if (length(x) == 0 || all(is.na(x))) NA_character_ else paste(x, collapse = ";"),
        character(1)
      )))
  }
  
  out
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

# ------------------
# 6. VOLCANO PLOTS 
# ------------------

contraste_meta <- list(
  HTNV_12h_vs_Control = list(titulo = "12h vs Control", sub = "Infección temprana (12h post-infección)",
                             col_up = "#F5A623", col_down = "#4C72B0"),
  HTNV_24h_vs_Control = list(titulo = "24h vs Control", sub = "Infección intermedia (24h post-infección)",
                             col_up = "#E8755A", col_down = "#4C72B0"),
  HTNV_36h_vs_Control = list(titulo = "36h vs Control", sub = "Infección tardía (36h post-infección)",
                             col_up = "#B0304A", col_down = "#4C72B0")
)

# N_LABELS_VOLCANO controla cuántos genes como máximo se etiquetan por volcano,
# por dirección (Up/Down) y por criterio (más significativos + más extremos
# en |log2FC|); tras deduplicar, el total visible suele quedar entre
# N_LABELS_VOLCANO y 2*N_LABELS_VOLCANO por dirección. Bájalo si sigue saturado.
N_LABELS_VOLCANO <- 12

plot_volcano <- function(res_df, titulo, sub, col_up, col_down,
                         highlight_genes = NULL, highlight_name = "Interés",
                         highlight_color = "#2CA02C",
                         n_labels = N_LABELS_VOLCANO) {
  
  df <- res_df %>% mutate(logp = -log10(padj + 1e-300))
  
  # Etiqueta la UNIÓN de: (a) los genes más significativos (menor padj) y
  # (b) los genes más extremos en magnitud (|log2FC| más alto). Solo por
  # significancia se quedaríam fuera genes muy desplazados en el eje X pero
  # con padj no-mínimo, que es justo lo que se ve como "extremos sin anotar".
  pick_top <- function(data, n) {
    if (nrow(data) == 0) return(data)
    bind_rows(
      data %>% slice_max(logp, n = n, with_ties = FALSE),
      data %>% slice_max(abs(log2FoldChange), n = n, with_ties = FALSE)
    ) %>% distinct(gene_id, .keep_all = TRUE)
  }
  
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
    colores <- setNames(c(highlight_color, col_up, col_down, "grey75"),
                        c(highlight_name, "Up", "Down", "NS"))
    top_lab <- bind_rows(
      pick_top(df %>% filter(grupo == highlight_name) %>%
                 filter(DEG_class != "NS" | abs(log2FoldChange) > 0.5), n_labels),
      pick_top(df %>% filter(DEG_class %in% c("Up", "Down")), n_labels)
    ) %>% distinct(gene_id, .keep_all = TRUE)
  } else {
    df <- df %>% mutate(grupo = DEG_class)
    colores <- c("Up" = col_up, "Down" = col_down, "NS" = "grey75")
    top_lab <- bind_rows(
      pick_top(df %>% filter(DEG_class == "Up"),   n_labels),
      pick_top(df %>% filter(DEG_class == "Down"), n_labels)
    ) %>% distinct(gene_id, .keep_all = TRUE)
  }
  
  ggplot(df, aes(x = log2FoldChange, y = logp, color = grupo)) +
    geom_point(data = df %>% filter(grupo == "NS"), alpha = 0.25, size = 1.3) +
    geom_point(data = df %>% filter(grupo != "NS"), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = colores, drop = FALSE, name = NULL) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(
      data = top_lab, aes(label = label_gene), size = 3.0,
      color = "black", fontface = "bold", max.overlaps = 20,
      force = 1.5, force_pull = 0.3, box.padding = 0.4,
      segment.size = 0.3, segment.color = "grey50", min.segment.length = 0.1,
      seed = 42
    ) +
    labs(title = titulo, subtitle = sub, x = "Log2 Fold Change",
         y = expression(-log[10](padj))) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", plot.subtitle = element_text(color = "grey40", size = 9))
}

# res_lista ya se construyó correctamente en la sección 5.
# Aquí se genera un volcano por contraste, guardando cada uno en
# plots_volcano.
plots_volcano <- list()
for (nombre in names(contrastes)) {
  m <- contraste_meta[[nombre]]
  plots_volcano[[nombre]] <- plot_volcano(
    res_lista[[nombre]],
    titulo = m$titulo, sub = m$sub, col_up = m$col_up, col_down = m$col_down
  )
}

p_panel <- (plots_volcano$HTNV_12h_vs_Control |
              plots_volcano$HTNV_24h_vs_Control |
              plots_volcano$HTNV_36h_vs_Control) +
  plot_annotation(
    title    = "Time-course infección HTNV — E-MTAB-11353 (mBMDM)",
    subtitle = sprintf("Cada tiempo post-infección vs control (mock, virus inactivado Co60, 0h) | |LFC| > 1, padj < 0.05 | etiquetas: top %d por padj",
                       N_LABELS_VOLCANO),
    theme    = theme(plot.title = element_text(size = 15, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

save_pub(paste0(output_dir, "Fig03_Volcano_Panel.png"), p_panel, width = 18, height = 6.5)

# -------------------------------
# 7. HEATMAP DEGs (3 CONTRASTES)
# -------------------------------
degs_union <- bind_rows(res_lista) %>%
  filter(DEG_class != "NS") %>%
  pull(gene_id) %>%
  unique()

cat(sprintf("DEGs unión de contrastes: %d genes\n", length(degs_union)))

if (length(degs_union) >= 5) {
  
  N_MAX <- 80
  if (length(degs_union) > N_MAX) {
    degs_union <- bind_rows(res_lista) %>%
      filter(gene_id %in% degs_union) %>%
      group_by(gene_id) %>%
      summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
      arrange(min_padj) %>%
      slice_head(n = N_MAX) %>%
      pull(gene_id)
    cat(sprintf("  Limitado a top %d por padj mínimo.\n", N_MAX))
  }
  
  gene_labels <- annotation %>%
    filter(gene_id %in% degs_union) %>%
    mutate(label = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id)) %>%
    dplyr::select(gene_id, label)
  
  mat_hm <- assay(vsd)[degs_union, ]
  rownames(mat_hm) <- make.unique(gene_labels$label[match(rownames(mat_hm), gene_labels$gene_id)])
  mat_scaled <- t(scale(t(mat_hm)))
  
  col_annot <- data.frame(Condition = coldata$condition, row.names = rownames(coldata))
  
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
  
  row_annot_ordered <- row_annot[degs_union, , drop = FALSE]
  rownames(row_annot_ordered) <- rownames(mat_scaled)
  
  deg_colors <- c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey90")
  ann_colors <- c(list(Condition = pal_cond),
                  setNames(rep(list(deg_colors), length(contrastes)), names(contrastes)))
  
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
    main     = "DEGs — E-MTAB-11353 (time-course HTNV, mBMDM)",
    filename = paste0(output_dir, "Fig04_heatmap_DEGs.png"),
    width = 13, height = max(10, n_genes * 0.18 + 3)
  )
}

# ---------------------------------------------------------
# 8. ENRIQUECIMIENTO FUNCIONAL (GO BP + KEGG, org.Mm.eg.db)
# ---------------------------------------------------------

run_enrichment <- function(res_df, nombre) {
  cat(sprintf("\nEnriquecimiento: %s\n", nombre))
  
  genes_up   <- res_df %>% filter(DEG_class == "Up",   !is.na(entrez)) %>% pull(entrez) %>% unique()
  genes_down <- res_df %>% filter(DEG_class == "Down", !is.na(entrez)) %>% pull(entrez) %>% unique()
  universe   <- res_df %>% filter(!is.na(entrez)) %>% pull(entrez) %>% unique()
  
  cat(sprintf("  Universo: %d | UP: %d | DOWN: %d\n", length(universe), length(genes_up), length(genes_down)))
  
  out <- list()
  
  if (length(genes_up) >= 5) {
    out$ego_up <- tryCatch(
      enrichGO(genes_up, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               readable = TRUE) %>% clusterProfiler::simplify(cutoff = 0.4, by = "p.adjust", select_fun = min),
      error = function(e) NULL
    )
    if (!is.null(out$ego_up) && nrow(as.data.frame(out$ego_up)) > 0)
      write.csv(as.data.frame(out$ego_up), paste0(results_dir, "GO_BP_up_", nombre, ".csv"), row.names = FALSE)
  }
  
  if (length(genes_down) >= 5) {
    out$ego_down <- tryCatch(
      enrichGO(genes_down, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               readable = TRUE) %>% clusterProfiler::simplify(cutoff = 0.4, by = "p.adjust", select_fun = min),
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

# ------------------------------------------------------------------------------
# 8a. COMPARACIÓN GO BP 
#     En vez de un dotplot por contraste (redundante entre sí y con muchos
#     términos que dicen casi lo mismo), se agrupan por similitud semántica
#     (rrvgo) TODOS los términos GO de los 3 contrastes juntos, se elige un
#     término representante por grupo, y se muestran solo los grupos no
#     redundantes en una única figura
# ------------------------------------------------------------------------------

N_GRUPOS_GO <- 10  # número máximo de grupos GO representantes a mostrar por figura

build_go_comparison <- function(direccion, threshold = 0.7, n_grupos_mostrar = N_GRUPOS_GO) {
  campo <- if (direccion == "up") "ego_up" else "ego_down"
  
  df_all <- bind_rows(lapply(names(contrastes), function(nombre) {
    e <- enrich_lista[[nombre]][[campo]]
    if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
    as.data.frame(e) %>% mutate(contraste = nombre)
  }))
  
  if (is.null(df_all) || length(unique(df_all$ID)) < 2) {
    cat(sprintf("  Sin suficientes términos GO %s para la comparación.\n", direccion))
    return(invisible(NULL))
  }
  
  ids_unicos <- unique(df_all$ID)
  simMatrix <- tryCatch(
    calculateSimMatrix(ids_unicos, orgdb = "org.Mm.eg.db", ont = "BP", method = "Rel"),
    error = function(e) { cat(sprintf("  rrvgo simMatrix falló (%s): %s\n", direccion, e$message)); NULL }
  )
  if (is.null(simMatrix) || nrow(simMatrix) < 2) return(invisible(NULL))
  
  # Score de cada término = su mejor (más significativo) padj visto en
  # cualquiera de los 3 contrastes, para que rrvgo elija bien el representante.
  scores_df <- df_all %>% group_by(ID) %>% summarise(score = max(-log10(p.adjust)), .groups = "drop")
  scores <- setNames(scores_df$score, scores_df$ID)
  scores <- scores[rownames(simMatrix)]
  
  reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = threshold, orgdb = "org.Mm.eg.db")
  n_terminos_orig <- nrow(reducedTerms)
  
  # Cada término original -> su grupo semántico (parentTerm representante)
  mapa_grupo <- reducedTerms %>% dplyr::select(go, parentTerm) %>% distinct()
  
  df_grouped <- df_all %>%
    left_join(mapa_grupo, by = c("ID" = "go")) %>%
    filter(!is.na(parentTerm))
  
  # Dentro de cada grupo semántico y cada contraste, solo se conserva el
  # término más significativo (evita puntos redundantes por celda).
  df_grouped <- df_grouped %>%
    group_by(parentTerm, contraste) %>%
    slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  n_grupos <- length(unique(df_grouped$parentTerm))
  cat(sprintf("  rrvgo GO_%s (comparación 3 tiempos): %d términos -> %d grupos no redundantes\n",
              direccion, n_terminos_orig, n_grupos))
  
  top_grupos <- df_grouped %>%
    group_by(parentTerm) %>%
    summarise(mejor_padj = min(p.adjust), .groups = "drop") %>%
    arrange(mejor_padj) %>%
    slice_head(n = n_grupos_mostrar) %>%
    pull(parentTerm)
  
  df_plot <- df_grouped %>%
    filter(parentTerm %in% top_grupos) %>%
    mutate(
      contraste_lbl = factor(contraste, levels = names(contrastes),
                             labels = sapply(contraste_meta[names(contrastes)], function(m) m$titulo)),
      logp = -log10(p.adjust)
    )
  
  # Ordena los términos en el eje Y por significancia global
  orden_terminos <- df_plot %>%
    group_by(parentTerm) %>%
    summarise(mejor_padj = min(p.adjust), .groups = "drop") %>%
    arrange(desc(mejor_padj)) %>%
    pull(parentTerm)
  df_plot$parentTerm <- factor(df_plot$parentTerm, levels = orden_terminos)
  
  titulo_dir <- if (direccion == "up") "UP" else "DOWN"
  
  p <- ggplot(df_plot, aes(x = contraste_lbl, y = parentTerm, size = Count, color = logp)) +
    geom_point() +
    scale_color_gradient(low = "#4C72B0", high = "#B0304A", name = expression(-log[10](padj))) +
    scale_size_continuous(name = "Genes", range = c(2, 9)) +
    labs(title = titulo_dir,
         subtitle = sprintf("%d términos GO originales agrupados en %d grupos no redundantes (rrvgo) | top %d mostrados",
                            n_terminos_orig, n_grupos, length(top_grupos)),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(color = "grey40", size = 9),
          axis.text.x = element_text(face = "bold"))
  
  write.csv(df_plot, paste0(results_dir, "GO_comparacion_", direccion, ".csv"), row.names = FALSE)
  
  # En vez de guardar la figura aquí, se devuelve el plot junto con los
  # metadatos necesarios para poder combinarlo con el otro sentido (up/down)
  # en una única figura conjunta.
  list(plot = p, df = df_plot, n_terminos_orig = n_terminos_orig,
       n_grupos = n_grupos, top_grupos = top_grupos)
}

go_comparacion_up   <- build_go_comparison("up")
go_comparacion_down <- build_go_comparison("down")

# --- Figura única combinando los paneles UP y DOWN (en vez de 2 figuras) ---
paneles_go <- list(UP = go_comparacion_up, DOWN = go_comparacion_down)
paneles_go <- paneles_go[!sapply(paneles_go, is.null)]

if (length(paneles_go) > 0) {
  p_go_combinado <- wrap_plots(lapply(paneles_go, `[[`, "plot"), nrow = 1) +
    plot_annotation(
      title = "GO BP — comparación time-course (E-MTAB-11353)",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  
  n_terminos_max <- max(sapply(paneles_go, function(x) length(x$top_grupos)))
  altura_fig <- max(5, 0.35 * n_terminos_max + 3)
  ancho_fig  <- 9 * length(paneles_go)
  
  save_pub(paste0(output_dir, "Fig06_GO_comparacion.png"), p_go_combinado,
           width = ancho_fig, height = altura_fig)
} else {
  cat("  Sin paneles GO (up/down) disponibles para generar la figura combinada.\n")
}

# ------------------------------------------------------------------------------
# 8b. GSEA [contraste principal (36h vs control (0h))]
# ------------------------------------------------------------------------------

gene_rank <- res_lista[[nombre_principal]] %>%
  filter(!is.na(entrez), !is.na(log2FoldChange), !is.na(pvalue)) %>%
  mutate(rank_score = sign(log2FoldChange) * (-log10(pvalue + 1e-300))) %>%
  arrange(desc(rank_score)) %>%
  dplyr::select(entrez, rank_score) %>%
  deframe()
gene_rank <- gene_rank[!duplicated(names(gene_rank))]

if (length(gene_rank) >= 100) {
  gsea_res <- tryCatch(
    gseGO(geneList = gene_rank, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
          ont = "BP", minGSSize = 15, maxGSSize = 500,
          pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE),
    error = function(e) { cat(sprintf("  GSEA falló: %s\n", e$message)); NULL }
  )
  
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    write.csv(as.data.frame(gsea_res),
              paste0(results_dir, "GSEA_GO_BP_", nombre_principal, ".csv"), row.names = FALSE)
    p_gsea <- dotplot(gsea_res, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(sprintf("GSEA GO BP — %s (E-MTAB-11353)", nombre_principal)) +
      theme_minimal(base_size = 10) +
      theme(axis.text.y = element_text(size = 8))
    save_pub(paste0(output_dir, "Fig05b_GSEA_", nombre_principal, ".png"), p_gsea, width = 14, height = 10)
  } else {
    cat("  GSEA: sin términos significativos (padj < 0.05).\n")
  }
} 

# ------------------
# 9. GUARDADO FINAL
# ------------------

write.csv(as.data.frame(counts(dds, normalized = TRUE)),
          paste0(results_dir, "normalized_counts.csv"))

resumen_degs <- lapply(names(contrastes), function(nombre) {
  res_lista[[nombre]] %>% count(DEG_class) %>% mutate(contraste = nombre)
}) %>% bind_rows() %>% pivot_wider(names_from = DEG_class, values_from = n, values_fill = 0)

write.csv(resumen_degs, paste0(results_dir, "Resumen_DEGs_por_contraste.csv"), row.names = FALSE)

cat("\n=== RESUMEN FINAL ===\n")
print(resumen_degs)
cat("\nFiguras en:", output_dir, "\n")
cat("Resultados en:", results_dir, "\n")

saveRDS(dds, "data/processed/dds_EMTAB11353.rds")
