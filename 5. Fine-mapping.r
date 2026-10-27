library(data.table)
library(susieR)
library(dplyr)
library(Rfast)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript run_susie_analysis.R <chromosome_number> <genotype_file> <phenotype_file> <output_file>", call. = FALSE)
}

chromosome_number <- as.numeric(args[1])
genotype_file <- args[2]
phenotype_file <- args[3]
output_file <- args[4]

### 1. Load genotype data
SNPINDELSVSTR <- fread(genotype_file)
SNPINDELSVSTR <- as.data.frame(SNPINDELSVSTR)
rownames(SNPINDELSVSTR) <- SNPINDELSVSTR$ID
SNPINDELSVSTR <- SNPINDELSVSTR[, -3:-9]  # Remove unnecessary columns
SNPINDELSVSTR <- na.omit(SNPINDELSVSTR)

### 2. Load phenotype data
BlockPhe <- fread(phenotype_file)
BlockPhe <- as.data.frame(BlockPhe)
rownames(BlockPhe) <- BlockPhe$pid

### 3. Initialize result storage
z3_ALLpip <- data.frame(
  position = numeric(0),
  PIP = numeric(0),
  CS_NumberID = numeric(0),
  PheCol4 = numeric(0),
  Var = character(0),
  stringsAsFactors = FALSE
)

# Filter phenotypes for target chromosome
BlockPhe <- BlockPhe[BlockPhe$`#Chr` == chromosome_number, ]

### 4. Main loop: iterate through each phenotype
for (NumberOfPhe in seq_len(nrow(BlockPhe))) {
  
  ## 4.1 Extract current phenotype
  Phe_test <- BlockPhe[NumberOfPhe, -1:-6]
  Phe_t <- as.data.frame(t(Phe_test))
  
  ## 4.2 Extract genotypes on current chromosome
  SNP_STR_SV_Selected <- SNPINDELSVSTR[SNPINDELSVSTR$`#CHROM` == chromosome_number, ]
  Pos <- as.numeric(BlockPhe[NumberOfPhe, 2])
  startPos <- Pos - 500000
  endPos <- Pos + 500000
  
  SNP_STR_SV_Selected_filtered <- SNP_STR_SV_Selected[
    SNP_STR_SV_Selected$POS >= startPos &
    SNP_STR_SV_Selected$POS <= endPos, -1:-2]
  
  ## 4.3 Transpose genotype data and check row names
  SNP_STR_SV_Selected_filtered_t <- as.data.frame(t(SNP_STR_SV_Selected_filtered))
  
  # Convert to numeric
  convert_to_numeric <- function(x) {
    if (is.numeric(x)) return(x)
    suppressWarnings({
      num_x <- as.numeric(as.character(x))
      if (all(is.na(num_x)) && !all(is.na(x))) {
        warning("Non-numeric values in column")
      }
      return(num_x)
    })
  }
  SNP_STR_SV_Selected_filtered_t[] <- lapply(SNP_STR_SV_Selected_filtered_t, convert_to_numeric)
  
  ## 4.4 Check if row names match
  if (!identical(rownames(SNP_STR_SV_Selected_filtered_t), rownames(Phe_t))) {
    stop("Row names do not match between genotype and phenotype data")
  }
  
  ## 4.5 Ensure consistent order
  SNP_STR_SV_Selected_filtered_t <- SNP_STR_SV_Selected_filtered_t[rownames(Phe_t), , drop = FALSE]
  
  ## 4.6 Run SuSiE analysis
  fitted <- susie(
    as.matrix(SNP_STR_SV_Selected_filtered_t),
    as.matrix(Phe_t),
    L = 10,
    verbose = TRUE
  )
  
  ## 4.7 Process and combine results
  length_of_cs <- length(fitted$sets$cs)
  
  for (CS_num in seq_len(length_of_cs)) {
    i <- fitted$sets$cs[[CS_num]]
    
    if (length(i) == 0) {
      cat("No SNPs in CS", CS_num, "\n")
      next
    }
    
    z3 <- cbind(i, fitted$pip[i])
    
    if (nrow(z3) > 0) {
      colnames(z3) <- c('position', 'PIP')
      z3_cur <- as.data.frame(z3)
      z3_cur$Var <- rownames(z3_cur)
      z3_cur$CS_NumberID <- CS_num
      z3_cur$PheCol4 <- BlockPhe[NumberOfPhe, 4]
      z3_ALLpip <- rbind(z3_ALLpip, z3_cur)
    } else {
      cat("No content in CS", CS_num, "\n")
    }
  }
}

### 5. Process final results
result_PIP <- z3_ALLpip %>%
  group_by(PheCol4) %>%
  ungroup()

### 6. Output results
fwrite(result_PIP, output_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, na = "NA")
cat("Analysis completed. Results saved to:", output_file, "\n")


