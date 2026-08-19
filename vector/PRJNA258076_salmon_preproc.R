# ==============================================================================
# PREPROCESAMIENTO PRJNA258076 — FASTQ → Counts (Salmon)
# ==============================================================================
# Este script:
# 1. Mapea archivos SRR.fastq a nombres de muestra (RR18, RR20, etc.)
# 2. Ejecuta Salmon para pseudoalineamiento rápido
# 3. Genera matriz de counts en formato esperado por DESeq2_PRJNA258076.R
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 0. CONFIGURACIÓN 
# ==============================================================================

# Mapeo SRR → muestra (ajusta si es diferente en tu descarga)
mapeo_srr <- data.frame(
  srr      = c("SRR1548195", "SRR1548706", "SRR1548707", "SRR1548708", "SRR1548709"),
  muestra  = c("RR18",       "RR20",       "RR29",       "RR30",       "RR31"),
  stringsAsFactors = FALSE
)

raw_fastq_dir  <- "data/raw/vector/PRJNA258076_RAW/"
salmon_index   <- "data/reference/Mus_musculus_ensembl_salmon_index/"
salmon_outdir  <- "data/processed/salmon_quant/"
counts_outdir  <- "data/raw/vector/PRJNA258076_RAW/"

dir.create(raw_fastq_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(salmon_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(counts_outdir, recursive = TRUE, showWarnings = FALSE)

# Verificar que salmon está instalado
if (system("which salmon > /dev/null 2>&1") != 0) {
  stop("salmon no está instalado. Instálalo con: conda install -c bioconda salmon")
}

# ==============================================================================
# 1. VALIDAR ARCHIVOS FASTQ
# ==============================================================================

cat("=== Validando archivos FASTQ ===\n")
archivos_faltantes <- c()

for (i in 1:nrow(mapeo_srr)) {
  srr <- mapeo_srr$srr[i]
  muestra <- mapeo_srr$muestra[i]
  
  # Buscar archivos: pueden ser .fastq, .fq, .fastq.gz, etc.
  fq_r1 <- list.files(raw_fastq_dir, 
                      pattern = sprintf("^%s.*_1\\.(fastq|fq)(\\.(gz|bz2|xz))?$", srr),
                      full.names = TRUE)
  fq_r2 <- list.files(raw_fastq_dir,
                      pattern = sprintf("^%s.*_2\\.(fastq|fq)(\\.(gz|bz2|xz))?$", srr),
                      full.names = TRUE)
}

if (length(archivos_faltantes) > 0) {
  cat("\n Archivos NO encontrados:\n")
  for (f in archivos_faltantes) cat(sprintf("  • %s\n", f))
}

# ==============================================================================
# 2. VALIDAR/DESCARGAR ÍNDICE SALMON
# ==============================================================================

if (!dir.exists(salmon_index)) {
  cat(sprintf("Índice no encontrado en %s\n", salmon_index))
  stop("Índice requerido")
}

# ==============================================================================
# 3. EJECUTAR SALMON PSEUDOALINEAMIENTO
# ==============================================================================


for (i in 1:nrow(mapeo_srr)) {
  srr <- mapeo_srr$srr[i]
  muestra <- mapeo_srr$muestra[i]
  
  fq_r1 <- list.files(
    raw_fastq_dir,
    pattern = sprintf("^%s.*_1\\.(fastq|fq)(\\.(gz|bz2|xz))?$", srr),
    full.names = TRUE
  )[1]
  
  fq_r2 <- list.files(
    raw_fastq_dir,
    pattern = sprintf("^%s.*_2\\.(fastq|fq)(\\.(gz|bz2|xz))?$", srr),
    full.names = TRUE
  )[1]
  
  outdir_muestra <- file.path(salmon_outdir, muestra)
  dir.create(outdir_muestra, recursive = TRUE, showWarnings = FALSE)
  
  # Si ya existe la cuantificación, no volver a ejecutarla
  quant_file <- file.path(outdir_muestra, "quant.sf")
  
  if (file.exists(quant_file)) {
    next
  }
  
  
  cmd <- sprintf(
    "salmon quant -i %s -l A -1 %s -2 %s -o %s --validateMappings --gcBias -p 4 2>/dev/null",
    salmon_index, fq_r1, fq_r2, outdir_muestra
  )
  
  ret <- system(cmd)
  
  if (ret == 0) {
    cat("✓\n")
  } else {
    warning(sprintf("Error ejecutando Salmon para %s", muestra))
  }
}
# ==============================================================================
# 4. CONVERTIR SALMON → FORMATO RSEM-compatible
# ==============================================================================

# Cargar anotación transcripto -> gen
tx2gene <- read.delim(
  "data/reference/tx2gene.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(mapeo_srr)) {
  
  muestra <- mapeo_srr$muestra[i]
  
  quant_file <- file.path(salmon_outdir, muestra, "quant.sf")
  
  if (!file.exists(quant_file)) {
    warning(sprintf("No existe %s", quant_file))
    next
  }
  
  df <- read.delim(quant_file)
  
  genes_results <-
    df %>%
    rename(
      transcript_id = Name,
      length = Length,
      effective_length = EffectiveLength,
      TPM = TPM,
      expected_count = NumReads
    ) %>%
    left_join(tx2gene, by = "transcript_id") %>%
    mutate(
      FPKM = NA_real_,
      IsoPct = NA_real_
    ) %>%
    select(
      gene_id,
      transcript_id,
      length,
      effective_length,
      expected_count,
      TPM,
      FPKM,
      IsoPct
    )
  
  outfile <- file.path(
    counts_outdir,
    paste0(muestra, ".genes.results")
  )
  
  write.table(
    genes_results,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

}

# ==============================================================================
# 5. VALIDACIÓN FINAL — intentar cargar con el script DESeq2
# ==============================================================================

muestras_test <- mapeo_srr$muestra
archivos_generados <- file.exists(file.path(counts_outdir, paste0(muestras_test, ".genes.results")))

if (all(archivos_generados)) {
  cat("Todos los archivos .genes.results generados correctamente\n")
} else {
  cat("Faltan archivos .genes.results para:\n")
  for (m in muestras_test[!archivos_generados]) {
    cat(sprintf("  • %s\n", m))
  }
}

