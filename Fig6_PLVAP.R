args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript generate_endothelial_metrics.R INPUT_RDATA OUTPUT_CSV", call. = FALSE)
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- args[[2]]

suppressPackageStartupMessages(library(Seurat))

plvap_state_genes <- c("PLVAP", "EMCN", "CA4", "RGCC", "ESM1", "ANGPT2", "APLN")
kdr_program_genes <- c("KDR", "ESM1", "ANGPT2", "APLN", "FLT1", "TEK")
vascular_state_genes <- unique(c(plvap_state_genes, kdr_program_genes))
endothelial_identity_genes <- c("PECAM1", "VWF", "CDH5", "ENG", "RAMP2")

zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

row_signature <- function(df, genes) {
  genes <- intersect(genes, colnames(df))
  if (!length(genes)) return(rep(NA_real_, nrow(df)))
  z <- as.data.frame(lapply(df[, genes, drop = FALSE], zscore))
  rowMeans(z, na.rm = TRUE)
}

loaded <- new.env(parent = emptyenv())
load(input_path, envir = loaded)
if (!exists("tiss_nonimmune", envir = loaded, inherits = FALSE)) {
  stop("The input RData does not contain 'tiss_nonimmune'.", call. = FALSE)
}

obj <- loaded$tiss_nonimmune
meta <- obj@meta.data
required_metadata <- c(
  "nonimmune_general_annotation", "analysis", "patient_id", "sample_type",
  "sample_name", "driver_gene"
)
missing_metadata <- setdiff(required_metadata, colnames(meta))
if (length(missing_metadata)) {
  stop(
    sprintf("Missing metadata columns: %s", paste(missing_metadata, collapse = ", ")),
    call. = FALSE
  )
}

endothelial_cells <- rownames(meta)[
  tolower(meta$nonimmune_general_annotation) == "endothelial"
]
meta_ec <- meta[endothelial_cells, , drop = FALSE]
endothelial_before_filter <- nrow(meta_ec)

stage_map <- c(naive = "TN", grouped_pr = "RD", grouped_pd = "PD")
meta_ec$stage <- factor(stage_map[meta_ec$analysis], levels = c("TN", "RD", "PD"))
meta_ec$patient <- sub("_NAT$", "", meta_ec$patient_id)
meta_ec$is_nat <- grepl("_NAT$", meta_ec$patient_id) |
  (!is.na(meta_ec$sample_type) & meta_ec$sample_type == "NAT")
meta_ec <- meta_ec[!meta_ec$is_nat & !is.na(meta_ec$stage), , drop = FALSE]
endothelial_cells <- rownames(meta_ec)

genes <- intersect(
  unique(c(vascular_state_genes, endothelial_identity_genes, "PLVAP")),
  rownames(obj)
)
expected_genes <- c(
  "PLVAP", "EMCN", "CA4", "RGCC", "ESM1", "ANGPT2", "APLN", "KDR",
  "FLT1", "TEK", "PECAM1", "VWF", "CDH5", "ENG", "RAMP2"
)
if (!identical(genes, expected_genes)) {
  stop("The source object's available ordered gene set differs from the verified input.", call. = FALSE)
}

expr <- as.matrix(
  GetAssayData(obj, assay = "RNA", slot = "data")[genes, endothelial_cells, drop = FALSE]
)

sample_keys <- unique(meta_ec[, c("sample_name", "patient", "stage", "driver_gene")])
sample_keys <- sample_keys[order(sample_keys$sample_name), ]
rows <- vector("list", nrow(sample_keys))
for (i in seq_len(nrow(sample_keys))) {
  key <- sample_keys[i, ]
  cells <- rownames(meta_ec)[meta_ec$sample_name == key$sample_name]
  idx <- match(cells, colnames(expr))
  values <- rowMeans(expr[, idx, drop = FALSE])
  detected <- rowMeans(expr[, idx, drop = FALSE] > 0)
  rows[[i]] <- cbind(
    key,
    n_ec = length(idx),
    as.data.frame(as.list(values)),
    plvap_positive_fraction = unname(detected["PLVAP"])
  )
}

sample_df <- do.call(rbind, rows)
sample_df$n_ec <- as.integer(sample_df$n_ec)
numeric_columns <- setdiff(
  colnames(sample_df),
  c("sample_name", "patient", "stage", "driver_gene")
)
sample_df[numeric_columns] <- lapply(sample_df[numeric_columns], as.numeric)
sample_df$plvap_state <- row_signature(sample_df, plvap_state_genes)
sample_df$kdr_program <- row_signature(sample_df, kdr_program_genes)
sample_df$vascular_state <- row_signature(sample_df, vascular_state_genes)

observed_stage_counts <- table(sample_df$stage)
if (
  nrow(meta) != 9830L ||
  endothelial_before_filter != 1098L ||
  nrow(meta_ec) != 857L ||
  sum(sample_df$n_ec) != 857L ||
  nrow(sample_df) != 30L ||
  ncol(sample_df) != 24L ||
  length(unique(sample_df$patient)) != 23L ||
  !identical(as.integer(observed_stage_counts), c(10L, 10L, 10L))
) {
  stop("Source-derived validation counts differ from the verified workflow.", call. = FALSE)
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(sample_df, output_path, row.names = FALSE)

cat(sprintf(
  "Generated %d sample rows from %d retained endothelial cells and %d patients.\n",
  nrow(sample_df), nrow(meta_ec), length(unique(sample_df$patient))
))

