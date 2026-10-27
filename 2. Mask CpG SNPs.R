#!/usr/bin/env Rscript

# Mask autosomal CpGs when either the C or adjacent G coordinate overlaps a SNP.
# Note that the reference genome coordinate systems for SNPs and CpGs must be consistent.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: mask_cpg_snps.R CpG_unmasked.bed.gz AllSNPList.txt CpG_masked.bed.gz")
}

raw_cpg       <- args[1]
snp_file      <- args[2]
output_cpg    <- args[3]

if (file.exists(output_cpg)) {
  stop("Refusing to overwrite existing output: ", output_cpg)
}

# -------------------------------------------------------------------------
# 1. Load SNP positions.
#    Expected format: chromosome:position, e.g. 1:123456
#    Chromosome prefixes such as "chr1" are normalized to "1".
# -------------------------------------------------------------------------

snp_lines <- readLines(snp_file, warn = FALSE)
snp_lines <- snp_lines[nzchar(trimws(snp_lines))]

snp_parts <- strsplit(snp_lines, ":", fixed = TRUE)
snp_chr   <- sub("^chr", "", vapply(snp_parts, `[`, character(1), 1),
                 ignore.case = TRUE)
snp_pos   <- as.integer(vapply(snp_parts, `[`, character(1), 2))

# Store SNP coordinates in a hash table for fast lookup.
snp_set <- new.env(hash = TRUE, parent = emptyenv())
for (key in paste0(snp_chr, ":", snp_pos)) {
  assign(key, TRUE, envir = snp_set)
}

# -------------------------------------------------------------------------
# 2. Read CpGs in chunks and mask SNP-overlapping sites.
#    The CpG coordinate represents the cytosine (C).
#    The adjacent guanine is therefore located at C + 1.
# -------------------------------------------------------------------------

in_con  <- gzfile(raw_cpg, "rt")
out_con <- gzfile(output_cpg, "wt")

on.exit({
  close(in_con)
  close(out_con)
}, add = TRUE)

autosomes <- as.character(1:29)

raw_count      <- 0L
masked_count   <- 0L
retained_count <- 0L

repeat {
  
  lines <- readLines(in_con, n = 100000L, warn = FALSE)
  if (!length(lines)) break
  
  fields <- strsplit(lines, "\t", fixed = TRUE)
  
  chr <- vapply(fields, `[`, character(1), 1)
  cpg <- as.integer(vapply(fields, `[`, character(1), 2))
  
  chr_key <- sub("^chr", "", chr, ignore.case = TRUE)
  is_autosomal <- chr_key %in% autosomes
  
  raw_count <- raw_count + sum(is_autosomal)
  
  # Mask a CpG if either its C or G coordinate overlaps a SNP.
  c_key <- paste0(chr_key, ":", cpg)
  g_key <- paste0(chr_key, ":", cpg + 1L)
  
  overlaps_snp <-
    vapply(c_key, exists, logical(1), envir = snp_set, inherits = FALSE) |
    vapply(g_key, exists, logical(1), envir = snp_set, inherits = FALSE)
  
  masked <- is_autosomal & overlaps_snp
  keep   <- is_autosomal & !overlaps_snp
  
  masked_count <- masked_count + sum(masked)
  
  # Renumber retained autosomal CpGs consecutively.
  n_keep <- sum(keep)
  
  if (n_keep > 0L) {
    ids <- retained_count + seq_len(n_keep)
    writeLines(
      paste(chr[keep], cpg[keep], ids, sep = "\t"),
      out_con
    )
    retained_count <- retained_count + n_keep
  }
}

# -------------------------------------------------------------------------
# 3. Write audit summary.
# -------------------------------------------------------------------------

audit_file <- paste0(output_cpg, ".audit.tsv")

writeLines(
  c(
    "category\tcount",
    paste0("raw_autosomal_cpgs\t", raw_count),
    paste0("masked_autosomal_cpgs\t", masked_count),
    paste0("retained_autosomal_cpgs\t", retained_count),
    paste0("source_raw_cpg\t", raw_cpg),
    paste0("source_snp_positions\t", snp_file),
    paste0("output_masked_cpg\t", output_cpg)
  ),
  audit_file
)

cat("raw_autosomal_cpgs\t",      raw_count,      "\n", sep = "")
cat("masked_autosomal_cpgs\t",   masked_count,   "\n", sep = "")
cat("retained_autosomal_cpgs\t", retained_count, "\n", sep = "")
