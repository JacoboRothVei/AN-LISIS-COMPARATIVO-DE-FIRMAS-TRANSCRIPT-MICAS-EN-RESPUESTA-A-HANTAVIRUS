#=============================================================
# FIRMA COMÚN HANTAVIRUS (BULK + scRNA) 
# Diagnóstico previo reveló:
#   - Severe y Control tienen muy pocos DEGs (73 y 117)
#   - Solapamiento máximo entre cualquier par = 446 (Acute∩scRNA)
#   - Severe ∩ Control = 0 genes → intersección estricta imposible
#   - scRNA tiene señal IFN mixta por tipo celular (esperado)
#   - Los 15 genes en n=3 son todos ISG/IFN canónicos (buena señal)
#=============================================================

R.version.string
rm(list=ls())
getwd()

#========================================
# BLOQUE 0: CONFIGURACIÓN Y LIBRERÍAS
#========================================
# Librerías
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(pheatmap)
library(RobustRankAggreg)
library(ggVennDiagram)
library(gridExtra)
library(patchwork)

# RobustRankAggreg 
rra_available <- requireNamespace("RobustRankAggreg", quietly = TRUE)


#Paths
gse133634       <- "results/human/GSE133634/DESeq2_results_full.csv"
gse158712_acute <- "results/human/GSE158712/DESeq2_acute_vs_recovered_full.csv"
gse158712_severe<- "results/human/GSE158712/DESeq2_severe_vs_moderate_full.csv"
gse161354       <- "results/human/GSE161354/DEGs_celulas_inmunes.csv"

output_dir  <- "figures/human/signature"
results_dir <- "results/human/signature"

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#Parámetros
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 0.5

# Anotacion para evitar conflictos
select <- dplyr::select
filter <- dplyr::filter

#========================================
# BLOQUE 1: CARGA DE DATOS
#========================================
leer_dataset <- function(ruta, nombre_estudio, scRNA = FALSE) {
  
  df <- read.csv(ruta, header = TRUE)
  
  if (scRNA) {
    df <- df %>%
      mutate(
        symbol         = trimws(as.character(gene)),
        log2FoldChange = avg_log2FC,
        padj           = p_val_adj
      )
  } else {
    df <- df %>%
      mutate(symbol = trimws(as.character(symbol)))
  }
  
  df %>%
    filter(
      !is.na(symbol), symbol != "",
      !is.na(padj),   padj < PADJ_CUTOFF,
      !is.na(log2FoldChange),
      abs(log2FoldChange) > LFC_CUTOFF
    ) %>%
    dplyr::select(symbol, log2FoldChange, padj) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    mutate(estudio = nombre_estudio)
}

bulk_acute  <- leer_dataset(gse158712_acute,  "Acute")
bulk_severe <- leer_dataset(gse158712_severe, "Severe")
bulk_ctrl   <- leer_dataset(gse133634,        "Control_vs_Infected")
sc_data     <- leer_dataset(gse161354,        "SingleCell", scRNA = TRUE)

cat(sprintf("  Acute: %d | Severe: %d | Control: %d | scRNA: %d\n",
            nrow(bulk_acute), nrow(bulk_severe),
            nrow(bulk_ctrl),  nrow(sc_data)))
#========================================
# BLOQUE 2: UNIVERSE LIMPIO
#========================================
universe_genes <- unique(c(
  bulk_acute$symbol, bulk_severe$symbol,
  bulk_ctrl$symbol,  sc_data$symbol
))
cat(sprintf("Universe (genes únicos medidos): %d\n\n", length(universe_genes)))

#========================================
# BLOQUE 3: FIRMA CORE — BULK CONSISTENCY
'''
un gen es "core hantavirus" si aparece como DEG
significativo en >= 2 de los 3 datasets bulk.
No exigimos same_direction porque bulk_ctrl puede tener
contraste invertido en algunos genes; usamos la mayoría.
scRNA se usa como evidencia de soporte, no como filtro.
'''
#========================================
master_df <- bind_rows(bulk_acute, bulk_severe, bulk_ctrl, sc_data)

bulk_df <- bind_rows(bulk_acute, bulk_severe, bulk_ctrl)

firma_core <- bulk_df %>%
  group_by(symbol) %>%
  summarise(
    n_bulk         = n_distinct(estudio),
    mean_lfc_bulk  = mean(log2FoldChange, na.rm = TRUE),
    # dirección mayoritaria (no unanimidad)
    direction      = ifelse(mean_lfc_bulk > 0, "UP", "DOWN"),
    n_up           = sum(log2FoldChange > 0),
    n_down         = sum(log2FoldChange < 0),
    min_padj_bulk  = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_bulk >= 2) %>%   # presente en >= 2 de 3 bulk
  arrange(min_padj_bulk)

cat(sprintf("Firma core bulk (>= 2/3 datasets): %d genes\n", nrow(firma_core)))

# Añadir soporte scRNA
sc_support <- sc_data %>%
  dplyr::select(symbol, lfc_sc = log2FoldChange, padj_sc = padj)

firma_core <- firma_core %>%
  left_join(sc_support, by = "symbol") %>%
  mutate(
    sc_support     = !is.na(lfc_sc),
    # ¿el scRNA va en la misma dirección que la mayoría bulk?
    sc_concordant  = sc_support & (sign(lfc_sc) == sign(mean_lfc_bulk)),
    # score de evidencia: combina magnitud, significancia y soporte sc
    evidence_score = abs(mean_lfc_bulk) *
      (-log10(min_padj_bulk)) *
      ifelse(sc_support, 1.25, 1)
  ) %>%
  arrange(desc(evidence_score))

write.csv(firma_core,
          file.path(results_dir, "firma_core_hantavirus.csv"),
          row.names = FALSE)

#========================================
# BLOQUE 4: RANK AGGREGATION (todos los datasets)
#========================================
score_ranked <- function(df) {
  df %>%
    mutate(score = abs(log2FoldChange) * (-log10(padj + 1e-300))) %>%
    arrange(desc(score)) %>%
    pull(symbol)
}

ranked_lists <- list(
  score_ranked(bulk_acute),
  score_ranked(bulk_severe),
  score_ranked(bulk_ctrl),
  score_ranked(sc_data)
)

if (rra_available) {
  rmat <- rankMatrix(ranked_lists)
  agg  <- aggregateRanks(rmat = rmat, method = "RRA")
  agg <- as.data.frame(agg)  
  # Corrección de Bonferroni sobre los p-valores RRA
  agg$Score_adj <- p.adjust(agg$Score, method = "BH")
  
  firma_rra <- agg %>%
    dplyr::filter(Score_adj < 0.05) %>%
    dplyr::rename(symbol = Name, pval_RRA = Score, padj_RRA = Score_adj)
  
  cat(sprintf("Genes firma RRA (BH < 0.05): %d\n", nrow(firma_rra)))
  
  write.csv(firma_rra,
            file.path(results_dir, "firma_RRA_hantavirus.csv"),
            row.names = FALSE)
  
  # Unión firma_core + RRA para enrichGO
  genes_enrichment <- union(firma_core$symbol, firma_rra$symbol)
} else {
  cat("RobustRankAggreg no disponible — usando sólo firma core\n")
  genes_enrichment <- firma_core$symbol
}

cat(sprintf("Genes para enrichment (core + RRA): %d\n\n", length(genes_enrichment)))

#========================================
# BLOQUE 5: enrichGO SOBRE FIRMA AMPLIADA
#========================================

ego <- enrichGO(
  gene          = genes_enrichment,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  universe      = universe_genes,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 5,
  maxGSSize     = 500
)

if (is.null(ego) || nrow(ego@result) == 0) {
  cat("AVISO: enrichGO sin resultados. Probando sin universe...\n")
  ego <- enrichGO(
    gene          = genes_enrichment,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    minGSSize     = 5
  )
}

ego_simp <- tryCatch(
  clusterProfiler::simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min),
  error = function(e) ego
)

n_sig <- nrow(ego_simp@result %>% filter(p.adjust < 0.05))
cat(sprintf("Términos GO-BP significativos: %d\n\n", n_sig))

write.csv(ego_simp@result,
          file.path(results_dir, "enrichGO_firma_ampliada.csv"),
          row.names = FALSE)

# Figura ORA
if (n_sig > 0) {
  p_ora <- dotplot(ego_simp, showCategory = 20, label_format = 50) +
    ggtitle("Firma hantavirus — ORA GO-BP") +
    theme_minimal(base_size = 11)
  ggsave(file.path(output_dir, "ORA_GO_v3.jpeg"),
         plot = p_ora, width = 10, height = 8, dpi = 300)
}

#========================================
# BLOQUE 6: GSEA INDEPENDIENTE POR DATASET
#========================================
convertir_ids <- function(df) {
  df %>%
    filter(!is.na(symbol), symbol != "") %>%
    distinct(symbol, .keep_all = TRUE) %>%
    { inner_join(.,
                 bitr(.$symbol, fromType="SYMBOL",
                      toType="ENTREZID", OrgDb=org.Hs.eg.db),
                 by = c("symbol"="SYMBOL")) }
}

run_gsea <- function(df, nombre, min_gs = 15) {
  
  df_m <- convertir_ids(df)
  if (nrow(df_m) == 0) { warning(paste("Sin genes:", nombre)); return(NULL) }
  
  gl <- setNames(df_m$log2FoldChange, df_m$ENTREZID)
  gl <- sort(gl, decreasing = TRUE)
  
  res <- tryCatch(
    gseGO(
      geneList      = gl,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      minGSSize     = min_gs,
      maxGSSize     = 500,
      eps           = 0,
      nPermSimple   = 10000
    ),
    error = function(e) { warning(paste("Error GSEA", nombre)); NULL }
  )
  
  if (is.null(res) || nrow(res@result) == 0) {
    cat(sprintf("  [%s] Sin resultados GSEA\n", nombre)); return(NULL)
  }
  
  n_sig <- nrow(res@result %>% filter(p.adjust < 0.05))
  cat(sprintf("  [%s] Pathways significativos: %d\n", nombre, n_sig))
  
  return(res)
}

#GSEA
gsea_acute  <- run_gsea(bulk_acute,  "Acute")
gsea_severe <- run_gsea(bulk_severe, "Severe")
gsea_ctrl   <- run_gsea(bulk_ctrl,   "Control_vs_Infected")
gsea_sc     <- run_gsea(sc_data,     "SingleCell", min_gs = 10)

# --------------------------------------------------------------------------
# Figura única comparativa GSEA — estilo lollipop divergente (igual que
# GO_comparativo.png): NES ya trae signo propio (activado/reprimido), así
# que no hace falta forzarlo como con GeneRatio en el ORA.
# --------------------------------------------------------------------------
N_TERMINOS_GSEA <- 10  # top N términos por dirección (activado/reprimido) y por dataset

clasificar_categoria_gsea <- function(term) {
  term_low <- tolower(term)
  case_when(
    str_detect(term_low, "interferon|ifn|antiviral|jak|stat|viral defense|defense response to virus") ~ "IFN / antiviral",
    str_detect(term_low, "cytokine|inflammat|neutrophil|monocyte|acute phase") ~ "Inflamación",
    str_detect(term_low, "endotheli|vascular|coagul|angiogenesis|platelet") ~ "Endotelio",
    str_detect(term_low, "nk cell|t cell|b cell|lymph|myeloid|innate|adaptive|leukocyte") ~ "Inmunidad celular",
    TRUE ~ "Otro"
  )
}

pal_gsea_categoria <- c(
  "IFN / antiviral"   = "#1B7FBF",
  "Inflamación"       = "#D95F02",
  "Endotelio"         = "#7570B3",
  "Inmunidad celular" = "#1B9E77",
  "Otro"              = "#999999"
)

extraer_df_gsea <- function(res_obj, nombre, n_top = N_TERMINOS_GSEA) {
  if (is.null(res_obj)) return(NULL)
  df <- res_obj@result
  if (nrow(df) == 0) return(NULL)
  
  df %>%
    filter(p.adjust < 0.05) %>%
    mutate(categoria = clasificar_categoria_gsea(Description),
           dataset   = nombre) %>%
    arrange(p.adjust) %>%
    group_by(sign(NES)) %>%
    slice_head(n = n_top) %>%
    ungroup()
}

df_gsea_todo <- bind_rows(
  extraer_df_gsea(gsea_acute,  "Acute"),
  extraer_df_gsea(gsea_severe, "Severe"),
  extraer_df_gsea(gsea_ctrl,   "Control_vs_Infected"),
  extraer_df_gsea(gsea_sc,     "SingleCell")
)

if (!is.null(df_gsea_todo) && nrow(df_gsea_todo) > 0) {
  # Orden compartido de términos en el eje Y: por categoría y luego por NES medio
  orden_terminos_gsea <- df_gsea_todo %>%
    group_by(Description, categoria) %>%
    summarise(orden_val = mean(NES), .groups = "drop") %>%
    arrange(categoria, orden_val) %>%
    pull(Description)
  
  df_gsea_todo <- df_gsea_todo %>%
    mutate(
      Description = factor(Description, levels = orden_terminos_gsea),
      dataset      = factor(dataset,
                            levels = c("Acute", "Severe", "Control_vs_Infected", "SingleCell"))
    )
  
  lim_x_gsea <- max(abs(df_gsea_todo$NES)) * 1.15
  
  p_gsea_comparativo <- ggplot(df_gsea_todo, aes(x = NES, y = Description)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey20") +
    geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description,
                     color = categoria), linewidth = 0.7) +
    geom_point(aes(size = setSize, color = categoria), alpha = 0.9) +
    scale_color_manual(values = pal_gsea_categoria, name = "Categoría funcional") +
    scale_size_continuous(name = "Set size", range = c(2, 8)) +
    scale_x_continuous(limits = c(-lim_x_gsea, lim_x_gsea)) +
    facet_wrap(~ dataset, nrow = 1) +
    labs(
      title    = "GSEA GO-BP — comparación entre datasets (firma hantavirus)",
      subtitle = sprintf("top %d términos por dirección (activado/reprimido) y dataset | padj<0.05",
                         N_TERMINOS_GSEA),
      x = "NES     \u2190  Reprimido            Activado  \u2192",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text      = element_text(face = "bold", size = 11),
      panel.spacing   = unit(1.2, "lines"),
      plot.subtitle   = element_text(color = "grey40", size = 9),
      legend.position = "right"
    )
  
  altura_fig_gsea <- max(6, 0.32 * length(orden_terminos_gsea) + 2.5)
  
  ggsave(file.path(output_dir, "GSEA_comparativo.jpeg"),
         plot = p_gsea_comparativo,
         width = 16, height = altura_fig_gsea, dpi = 300, limitsize = FALSE)
} else {
  cat("  Sin términos GSEA disponibles para generar la figura comparativa.\n")
}

#========================================
# BLOQUE 7: PATHWAYS COMUNES ENTRE DATASETS
#========================================
get_sig_terms <- function(obj, p_cut = 0.05) {
  if (is.null(obj)) return(character(0))
  obj@result %>% filter(p.adjust < p_cut) %>% pull(Description)
}

terms_list <- list(
  Acute   = get_sig_terms(gsea_acute),
  Severe  = get_sig_terms(gsea_severe),
  Control = get_sig_terms(gsea_ctrl),
  scRNA   = get_sig_terms(gsea_sc)
)

all_terms  <- unlist(terms_list)
term_count <- table(all_terms)

common_2 <- names(term_count[term_count >= 2])
common_3 <- names(term_count[term_count >= 3])

cat(sprintf("\nPathways GSEA comunes (>=2 datasets): %d\n", length(common_2)))
cat(sprintf("Pathways GSEA comunes (>=3 datasets): %d\n", length(common_3)))

write.csv(
  data.frame(
    pathway    = names(term_count),
    n_datasets = as.integer(term_count)
  ) %>% arrange(desc(n_datasets)),
  file.path(results_dir, "pathways_comunes_GSEA_v3.csv"),
  row.names = FALSE
)

#========================================
# BLOQUE 8: FIGURA PRINCIPAL — FIRMA BIOLÓGICA
#========================================
keywords  <- c("interferon","IFN","antiviral","JAK","STAT",
               "cytokine","inflammat","neutrophil","monocyte",
               "endotheli","vascular","NK","T cell","B cell",
               "lymph","myeloid","innate","defense","viral")
kw_pat <- paste(keywords, collapse = "|")

# Combinar términos ORA + GSEA relevantes
gsea_all <- bind_rows(
  if (!is.null(gsea_acute))  gsea_acute@result  %>% mutate(dataset="Acute"),
  if (!is.null(gsea_severe)) gsea_severe@result %>% mutate(dataset="Severe"),
  if (!is.null(gsea_ctrl))   gsea_ctrl@result   %>% mutate(dataset="Control"),
  if (!is.null(gsea_sc))     gsea_sc@result     %>% mutate(dataset="scRNA")
)

firma_bio <- bind_rows(
  # ORA
  ego_simp@result %>%
    filter(p.adjust < 0.05,
           grepl(kw_pat, Description, ignore.case = TRUE)) %>%
    transmute(Description, pvalue = p.adjust, NES = NA_real_,
              source = "ORA", dataset = "all"),
  # GSEA
  if (nrow(gsea_all) > 0)
    gsea_all %>%
    filter(p.adjust < 0.05,
           grepl(kw_pat, Description, ignore.case = TRUE)) %>%
    transmute(Description, pvalue = p.adjust, NES,
              source = "GSEA", dataset)
) %>%
  group_by(Description) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(pvalue)

cat(sprintf("\nTérminos firma biológica integrada: %d\n", nrow(firma_bio)))

write.csv(firma_bio,
          file.path(results_dir, "firma_biologica_v3.csv"),
          row.names = FALSE)

if (nrow(firma_bio) > 0) {
  
  top_bio <- firma_bio %>%
    slice_head(n = 20) %>%
    mutate(
      logP      = -log10(pvalue),
      categoria = case_when(
        grepl("interferon|IFN|antiviral|JAK|STAT|viral defense",
              Description, ignore.case = TRUE)          ~ "IFN / antiviral",
        grepl("cytokine|inflammat|neutrophil|monocyte",
              Description, ignore.case = TRUE)          ~ "Inflamación",
        grepl("endotheli|vascular|coagul",
              Description, ignore.case = TRUE)          ~ "Endotelio",
        grepl("NK|T cell|B cell|lymph|myeloid|innate|adaptive",
              Description, ignore.case = TRUE)          ~ "Inmunidad celular",
        TRUE ~ "Otro"
      )
    )
  
  colores <- c("IFN / antiviral"   = "#1B7FBF",
               "Inflamación"       = "#D95F02",
               "Endotelio"         = "#7570B3",
               "Inmunidad celular" = "#1B9E77",
               "Otro"              = "#999999")
  
  p_final <- ggplot(top_bio,
                    aes(x = reorder(Description, logP), y = logP,
                        fill = categoria, alpha = source)) +
    geom_col() +
    scale_alpha_manual(values = c("ORA"=1, "GSEA"=0.75), name="Método") +
    scale_fill_manual(values = colores, name = "Proceso") +
    coord_flip() +
    labs(
      title    = "Hantavirus response signature v3",
      subtitle = sprintf("Genes core: %d | Enrichment genes: %d",
                         nrow(firma_core), length(genes_enrichment)),
      x = NULL, y = expression(-log[10](p.adjust))
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title   = element_text(face="bold", size=13),
          axis.text.y  = element_text(size=9),
          legend.position = "right")
  
  ggsave(file.path(output_dir, "Figura_firma_v3.jpeg"),
         plot = p_final, width = 13, height = 7, dpi = 300)
}

#========================================
# BLOQUE 9: HEATMAP FIRMA CORE
# LFC de los top genes core en cada dataset
#========================================
top_core_genes <- firma_core %>%
  slice_head(n = 30) %>%
  pull(symbol)

heatmap_df <- master_df %>%
  dplyr::filter(symbol %in% top_core_genes) %>%
  dplyr::select(symbol, estudio, log2FoldChange) %>%
  pivot_wider(names_from = estudio, values_from = log2FoldChange)

# Guardar para inspección
write.csv(heatmap_df,
          file.path(results_dir, "heatmap_firma_core.csv"),
          row.names = FALSE)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  mat <- heatmap_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()
  mat_viz <- mat; mat_viz[is.na(mat_viz)] <- 0 # Reemplazar NA por 0 sólo para visualización
  
  pheatmap(
    mat_viz,
    color            = colorRampPalette(c("#2166AC","white","#D6604D"))(100),
    cluster_cols     = FALSE,
    fontsize_row     = 8,
    main             = "LFC por dataset — firma core hantavirus\n(NA = no DEG en ese estudio)",
    filename         = file.path(output_dir, "Heatmap_firma_core.jpeg"),
    width = 7, height = 9
  )
}


#========================================
# BLOQUE 10: DIAGRAMAS DE VENN
#========================================


# Receptores y genes de interés Hantavirus
genes_hanta_receptores <- c(
  "CD55", "PCDH1", "CGAS", "TLR3", "TLR8",
  "STING1", "FLT1", "CDKN1A", "TLR4", "CLEC4M",
  "ITGAV", "ITGB3", "DAG1"
)

# Listas de genes
top100_lists <- list(
  "Acute\n(GSE158712)"  = bulk_acute  %>% arrange(padj) %>% slice_head(n = 100) %>% pull(symbol),
  "Severe\n(GSE158712)" = bulk_severe %>% arrange(padj) %>% slice_head(n = 100) %>% pull(symbol),
  "HTNV\n(GSE133634)"   = bulk_ctrl   %>% arrange(padj) %>% slice_head(n = 100) %>% pull(symbol),
  "scRNA\n(GSE161354)"  = sc_data     %>% arrange(padj) %>% slice_head(n = 100) %>% pull(symbol)
)

recept_lists <- list(
  "Acute\n(GSE158712)"  = bulk_acute  %>% dplyr::filter(symbol %in% genes_hanta_receptores) %>% pull(symbol),
  "Severe\n(GSE158712)" = bulk_severe %>% dplyr::filter(symbol %in% genes_hanta_receptores) %>% pull(symbol),
  "HTNV\n(GSE133634)"   = bulk_ctrl   %>% dplyr::filter(symbol %in% genes_hanta_receptores) %>% pull(symbol),
  "scRNA\n(GSE161354)"  = sc_data     %>% dplyr::filter(symbol %in% genes_hanta_receptores) %>% pull(symbol)
)

for (nm in names(recept_lists)) {
  cat(sprintf("  %s: %s\n", nm,
              ifelse(length(recept_lists[[nm]]) == 0, "ninguno",
                     paste(recept_lists[[nm]], collapse = ", "))))
}

# Función principal 
venn_con_tabla <- function(gene_lists, titulo, subtitulo, color_high, ruta) {
  
  # Venn plot
  p_venn <- ggVennDiagram(
    gene_lists, label = "count", label_alpha = 0, edge_size = 0.5
  ) +
    scale_fill_gradient(low = "#FFF5EB", high = color_high) +
    scale_color_manual(values = rep("grey40", length(gene_lists))) +
    labs(title = titulo, subtitle = subtitulo, fill = "N genes") +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey40", size = 7))
  
  # Extraer solapamientos exclusivos (SOLO solapamientos entre 2 o más conjuntos)
  
  sets <- gene_lists
  nms  <- names(sets)
  
  overlap_df <- bind_rows(lapply(2:length(nms), function(i) {  # empieza en 2
    combs <- combn(nms, i, simplify = FALSE)
    bind_rows(lapply(combs, function(grp) {
      genes_comunes <- Reduce(intersect, sets[grp])
      no_grp <- setdiff(nms, grp)
      if (length(no_grp) > 0) {
        genes_excluir <- unique(unlist(sets[no_grp]))
        genes_comunes <- setdiff(genes_comunes, genes_excluir)
      }
      if (length(genes_comunes) > 0) {
        data.frame(
          Conjuntos = paste(grp, collapse = " ∩ "),
          N         = length(genes_comunes),
          Genes     = paste(sort(genes_comunes), collapse = ", ")
        )
      }
    }))
  })) %>% dplyr::filter(N > 0) %>% arrange(desc(N))
  
  # Guardar CSV
  write.csv(overlap_df,
            file.path(results_dir, basename(gsub(".jpeg", "_tabla_genes.csv", ruta))),
            row.names = FALSE)
  print(overlap_df)
  
  # Tabla visual + Venn combinados
  if (nrow(overlap_df) > 0) {
    overlap_df_viz <- overlap_df %>%
      mutate(
        Conjuntos = str_wrap(Conjuntos, width = 30),
        Genes     = str_wrap(Genes,     width = 60)
      )
    
    tabla <- tableGrob(overlap_df_viz, rows = NULL,
                       theme = ttheme_minimal(base_size = 8))
    
    p_combined <- gridExtra::arrangeGrob(
      p_venn, tabla,
      ncol    = 1,
      heights = c(3, nrow(overlap_df_viz) * 0.4 + 0.5)
    )
    
    ggsave(ruta, plot = p_combined,
           width  = 12,
           height = 8 + nrow(overlap_df_viz) * 0.5,
           dpi    = 300)
  } else {
    ggsave(ruta, plot = p_venn, width = 10, height = 8, dpi = 300)
  }
  
  return(invisible(p_venn))
}

# Llamadas
venn_con_tabla(
  gene_lists = top100_lists,
  titulo     = "Solapamiento Top 100 DEGs por dataset",
  subtitulo  = "Genes ordenados por padj, top 100 de cada estudio",
  color_high = "#2171B5",
  ruta       = file.path(output_dir, "Venn_top100_datasets.jpeg")
)

venn_con_tabla(
  gene_lists = recept_lists,
  titulo     = "Receptores Hantavirus detectados como DEG",
  subtitulo  = paste("Genes evaluados:", paste(genes_hanta_receptores, collapse = ", ")),
  color_high = "#D94801",
  ruta       = file.path(output_dir, "Venn_receptores_hantavirus.jpeg")
)

