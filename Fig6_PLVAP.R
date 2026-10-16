#!/usr/bin/env Rscript
# Usage: Rscript Fig6_PLVAP.R input.rds output_directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L)
  stop('Usage: Rscript Fig6_PLVAP.R input.rds [output_directory]')
suppressPackageStartupMessages(library(SeuratObject))
input <- normalizePath(args[1], mustWork = TRUE)
out <- if (length(args) == 2L) args[2] else 'PLVAP_results'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# RData support allows the same workflow to use an archived Seurat workspace.
if (grepl('\\.(RData|rda)$', input, ignore.case = TRUE)) {
  env <- new.env(); load(input, envir = env)
  candidates <- Filter(function(n) inherits(env[[n]], 'Seurat'), ls(env))
  if ('tiss_nonimmune' %in% candidates) obj <- env$tiss_nonimmune
  else if (length(candidates) == 1L) obj <- env[[candidates]]
  else stop('RData must contain tiss_nonimmune or one unambiguous Seurat object.')
} else obj <- readRDS(input)
if (!inherits(obj, 'Seurat')) stop('Input must be an annotated Seurat object.')
meta <- obj[[]]
required <- c('patient_id', 'sample_name', 'sample_type', 'analysis')
if (!all(required %in% names(meta)))
  stop('Missing metadata: ', paste(setdiff(required, names(meta)), collapse = ', '))
# The input must be the nonimmune subset used for this analysis.
if ('immune_annotation' %in% names(meta)) {
  labels <- unique(tolower(trimws(as.character(meta$immune_annotation))))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (any(!labels %in% c('non-immune', 'nonimmune', 'non_immune', 'non immune')))
    stop('Input contains cells outside the nonimmune subset. Supply the nonimmune Seurat object.')
}
# Read the existing normalized RNA data without renormalizing or using scale.data.
if (packageVersion('SeuratObject') >= '5.0.0') {
  layers <- Layers(obj[['RNA']])
  if (!'data' %in% layers)
    stop('RNA assay needs a single normalized data layer. Join split data layers before running.')
  expression <- LayerData(obj, assay = 'RNA', layer = 'data')
} else expression <- GetAssayData(obj, assay = 'RNA', slot = 'data')
genes <- c('PLVAP', 'PECAM1', 'PTPRC')
if (!all(genes %in% rownames(expression))) stop('RNA data must contain PLVAP, PECAM1, PTPRC.')
if (!all(rownames(meta) %in% colnames(expression))) stop('Metadata/expression cell IDs do not match.')
v <- as.matrix(expression[genes, rownames(meta), drop = FALSE])
if (any(!is.finite(v)) || any(v < 0)) stop('Expected finite nonnegative RNA log-normalized expression.')
stage_map <- c(naive = 'TN', grouped_pr = 'RD', grouped_pd = 'PD')
patient <- as.character(meta$patient_id)
stage <- unname(stage_map[as.character(meta$analysis)])
keep <- !is.na(stage) & !is.na(patient) & nzchar(patient) &
  !grepl('_NAT$', patient) & (is.na(meta$sample_type) | meta$sample_type != 'NAT')
if (any(is.na(meta$sample_name[keep]) | !nzchar(as.character(meta$sample_name[keep]))))
  stop('Retained cells have missing biopsy/sample_name values.')
a <- data.frame(cell_barcode = rownames(meta)[keep],
  patient = sub('_NAT$', '', patient[keep]), stage = stage[keep],
  biopsy = as.character(meta$sample_name[keep]),
  PLVAP = v['PLVAP', keep], PECAM1 = v['PECAM1', keep], PTPRC = v['PTPRC', keep])
a$denominator <- as.integer(a$PECAM1 > 2 & a$PTPRC == 0)
a$numerator <- as.integer(a$denominator == 1L & a$PLVAP > 2)
b <- aggregate(cbind(numerator, denominator) ~ patient + stage + biopsy, a, sum)
b$fraction <- ifelse(b$denominator > 0, b$numerator / b$denominator, NA_real_)
valid <- b[!is.na(b$fraction), ]
if (!nrow(valid)) stop('No biopsies have eligible cells.')
p <- aggregate(fraction ~ patient + stage, valid, mean)
p <- p[order(p$patient, p$stage), ]
ns <- sapply(c('TN', 'RD', 'PD'), function(st) sum(p$stage == st))
if (any(ns < 2)) stop('At least two patient-stage values per group are required.')
# Patient-aligned version preserves repeated patient identities and missing stages.
aligned <- reshape(p, idvar = 'patient', timevar = 'stage', direction = 'wide')
names(aligned) <- sub('^fraction\\.', '', names(aligned))
aligned <- aligned[order(aligned$patient), c('patient', 'TN', 'RD', 'PD')]
# Compact three-column version: each column lists its own patients alphabetically.
# Rows across columns do not represent matched patients.
compact <- as.data.frame(setNames(lapply(c('TN', 'RD', 'PD'), function(st) {
  z <- p$fraction[p$stage == st]; c(z, rep(NA_real_, max(ns) - length(z)))
}), c('TN', 'RD', 'PD')))
tests <- do.call(rbind, lapply(c('TN', 'PD'), function(st) {
  x <- p$fraction[p$stage == st]; y <- p$fraction[p$stage == 'RD']
  test <- t.test(x, y, paired = FALSE, var.equal = TRUE, alternative = 'two.sided')
  data.frame(comparison = paste(st, 'RD', sep = '-'), n_group = length(x), n_RD = length(y),
    mean_group = mean(x), mean_RD = mean(y), difference = mean(x) - mean(y),
    t = unname(test$statistic), df = unname(test$parameter), raw_p = test$p.value,
    difference_CI95_lower = test$conf.int[1], difference_CI95_upper = test$conf.int[2])
}))
options(digits = 17)
write.csv(compact, file.path(out, 'TN_RD_PD_patient_fractions.csv'), row.names = FALSE, na = '')
write.csv(aligned, file.path(out, 'patient_fractions_with_IDs.csv'), row.names = FALSE, na = '')
write.csv(b, file.path(out, 'biopsy_counts_and_fractions.csv'), row.names = FALSE, na = '')
write.csv(tests, file.path(out, 'unpaired_Student_t_tests.csv'), row.names = FALSE)
write.csv(a, file.path(out, 'cell_selection_audit.csv'), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'session_info.txt'))
print(tests)
cat('\nSaved results to:', normalizePath(out), '\n')
