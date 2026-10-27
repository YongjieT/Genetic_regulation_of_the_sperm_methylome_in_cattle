#!/usr/bin/env Rscript
# meQTL–GWAS colocalization using coloc.abf
# PP.H4 > 0.60

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(doParallel)
  library(foreach)
})

# Usage: Rscript "10. Colocalization.R" trait cores pairs.tsv gwas.tsv meQTL_ids.txt GWAS_N [output_dir]

# 1. Parameters
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 6:7)

trait <- args[1]
cores <- as.integer(args[2])
N_gwas <- as.numeric(args[6])

script_dir <- dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
))
out <- if (length(args) == 7) args[7] else
  file.path(script_dir, "4_analysis", "coloc_primary_runs", trait)

if (dir.exists(out)) stop("Output directory already exists: ", out)
dir.create(out, recursive = TRUE)

# 2. Read and merge data
pairs <- fread(args[3])
gwas <- fread(args[4])
ids <- as.character(fread(args[5], header = FALSE)[[1]])

if (!length(ids) || anyNA(ids) || anyDuplicated(ids))
  stop("Missing or duplicate block IDs")

# Support prepared and legacy column names
if ("chr_pos" %in% names(pairs)) pairs[, variant_id := chr_pos]
if ("chr_pos" %in% names(gwas)) gwas[, chr_pos_ARSUCD1.2 := chr_pos]
if ("pval_gwas" %in% names(gwas)) gwas[, Pval_GWAS := pval_gwas]

# Retain target blocks and all shared SNPs
pairs <- pairs[phenotype_id %in% ids]
merged <- merge(pairs, gwas,
                by.x = "variant_id", by.y = "chr_pos_ARSUCD1.2",
                sort = FALSE)
grouped <- split(merged, by = "phenotype_id")

# 3. Parallel coloc.abf
cl <- makeCluster(cores)
registerDoParallel(cl)

results <- foreach(id = ids, .packages = c("coloc", "data.table")) %dopar% {
  d <- grouped[[id]]
  n <- if (is.null(d)) 0L else nrow(d)

  row <- data.table(
    phenotype_id = id, trait = trait,
    status = "NO_OVERLAP", nsnps = n,
    PP_H0 = NA_real_, PP_H1 = NA_real_,
    PP_H2 = NA_real_, PP_H3 = NA_real_, PP_H4 = NA_real_,
    max_SNP_PP_H4 = NA_real_, max_SNP_chr_pos = NA_character_
  )

  if (n >= 2L) {
    fit <- tryCatch(
      coloc.abf(
        dataset1 = list(snp = d$variant_id, pvalues = d$pval_nominal,
                        type = "quant", N = 125),
        dataset2 = list(snp = d$variant_id, pvalues = d$Pval_GWAS,
                        type = "quant", N = N_gwas),
        MAF = d$MAF
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      row[, status := "COLOC_ERROR"]
    } else {
      lead <- which.max(fit$results$SNP.PP.H4)

      row[, (paste0("PP_H", 0:4)) :=
            as.list(as.numeric(fit$summary[paste0("PP.H", 0:4, ".abf")]))]

      row[, `:=`(
        status = "OK",
        max_SNP_PP_H4 = fit$results$SNP.PP.H4[lead],
        max_SNP_chr_pos = fit$results$snp[lead]
      )]

      fwrite(fit$results,
             file.path(out, paste0("colocalization_", id, "_", trait, ".tsv")),
             sep = "\t")
    }
  }
  row
}

stopCluster(cl)

# 4. Save colocalization results
summary <- rbindlist(results)
fwrite(summary, file.path(out, "colocalization_summary.tsv"), sep = "\t")

if (any(summary$status == "COLOC_ERROR"))
  stop("Some coloc analyses failed; inspect colocalization_summary.tsv")

# 5. Select colocalized blocks (PP.H4 > 0.60)
selected <- summary[status == "OK" & !is.na(PP_H4) & PP_H4 > 0.60]
fwrite(selected, file.path(out, "coloc_regions_PP4_gt_0p6.tsv"), sep = "\t")

# 6. Deduplicate representative SNPs
representatives <- selected[, .(
  associated_blocks = paste(phenotype_id, collapse = ";"),
  number_of_blocks = .N
), by = .(trait, max_SNP_chr_pos)]

fwrite(representatives,
       file.path(out, "representative_snps_deduplicated.tsv"),
       sep = "\t")

# 7. Trait-level summary
fwrite(data.table(
  trait = trait,
  requested_blocks = length(ids),
  evaluable_blocks = sum(summary$status == "OK"),
  selected_blocks = nrow(selected),
  representative_snps = nrow(representatives)
), file.path(out, "trait_summary.tsv"), sep = "\t")

writeLines("COMPLETE", file.path(out, "COMPLETE"))
