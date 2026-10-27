#!/usr/bin/env bash
#semen quality as an example
set -euo pipefail

# 1. Paths and parameters
[[ $# -le 2 ]] || { echo "Usage: bash $0 [trait] [output.csv]" >&2; exit 2; }

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project="${MWAS_PROJECT:-$script_dir/4_analysis}"
trait="${1:-VE}"
python="$(command -v "${PYTHON_BIN:-python3}")"

case "$trait" in
    VE|MS|NMSP|NSP|SC) ;;
    *) echo "Unsupported trait: $trait" >&2; exit 2 ;;
esac

model_dir="$project/models/full_cis"
model_db="$model_dir/model.db"
covariance="$model_dir/covariance.txt.gz"
gwas="$project/gwas/SpermQuality/$trait/gwas.tsv.gz"
output="${2:-$project/4_analysis/direct_strategy_B/${trait}_MWAS_results.csv}"

# 2. Check inputs and output
[[ ! -e "$output" ]] || { echo "Output exists: $output" >&2; exit 1; }

for file in "$python" "$software" "$model_db" "$covariance" "$gwas"; do
    [[ -f "$file" && -r "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
mkdir -p "$(dirname "$output")"

# 3. SPrediXcan analysis
# Model training: 125 cattle
# Covariates: GPC1-5 and DNAmPC1-8 regressed out before CV.
# cis-window: meBlock start-500 kb to end+500 kb; original genotype filtering.
# GWAS: BETA for A1; A2 is the alternative allele; use BETA/SE.
# Downstream: restrict to 29,581 fine-mapped targets (target_universe.tsv).

"$python" SPrediXcan.py \
    --model_db_path "$model_db" \
    --covariance "$covariance" \
    --gwas_file "$gwas" \
    --snp_column SNP \
    --effect_allele_column A1 --non_effect_allele_column A2 \
    --beta_column BETA --se_column SE \
    --keep_non_rsid --stream_covariance --additional_output --throw \
    --output_file "$output"
