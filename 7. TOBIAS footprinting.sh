#!/usr/bin/env bash

# Set PROJECT_ROOT and CONDA_SH; configure the Slurm partition/resources.
# Six independent samples: ATACorrect -> ScoreBigwig -> BINDetect.

set -euo pipefail
export LC_ALL=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg

# Activate the existing TOBIAS environment.
source "${CONDA_SH:Set CONDA_SH to the Conda initialization script}"
conda activate "${TOBIAS_ENV:-tobias}"

# Require the TOBIAS version used in the original analysis.
version=$(TOBIAS --version 2>&1)
echo "$version"
if [[ ! "$version" =~ (^|[^0-9])0\.17\.5([^0-9]|$) ]]; then
    echo "ERROR: Expected TOBIAS 0.17.5" >&2
    exit 1
fi

# Read-only inputs: BAM, ARS-UCD1.2 FASTA, peak union and JASPAR2024 motifs.
project="${PROJECT_ROOT}"
reference="$project/01_reference"
sample="bull_${SLURM_ARRAY_TASK_ID}"
cores="${SLURM_CPUS_PER_TASK}"

bam="$project/$sample/usable.bam"
genome="$reference/genome.fa"
peaks="$reference/union_peaks.bed"
motifs="$reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

# Use one new run/sample directory per submission.
run="$project/tobias_slurm_${SLURM_ARRAY_JOB_ID}"
output="$run/$sample"
stage="input validation"
trap 'echo "FAILED STEP: $stage; SAMPLE: $sample; OUTPUT: $output; inspect Slurm error log before rerunning." >&2' ERR

for input in "$bam" "$genome" "$genome.fai" "$peaks" "$motifs"; do
    if [[ ! -s "$input" ]]; then
        echo "ERROR: Missing or empty input: $input" >&2
        exit 1
    fi
done
if [[ ! -s "$bam.bai" && ! -s "${bam%.bam}.bai" && ! -s "$bam.csi" ]]; then
    echo "ERROR: Missing BAM index: $bam" >&2
    exit 1
fi
mkdir -p "$run"
# Refuse to overwrite an existing sample directory for this submission.
mkdir "$output"
mkdir "$output/01_atacorrect" "$output/02_footprint" "$output/03_bindetect"
printf '%s\n' "$version" > "$output/software_version.txt"
printf 'sample\tbam\tgenome\tpeaks\tmotifs\n%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$bam" "$genome" "$peaks" "$motifs" > "$output/inputs.tsv"

# Trace executed commands in the Slurm log.
set -x

# 1. Correct Tn5 sequence bias.
# usable.bam must be unshifted; ATACorrect applies the +4/-5 shift once.
# Shared peak union; no blacklist; remaining options use TOBIAS defaults.
stage="ATACorrect"
TOBIAS ATACorrect \
    --bam "$bam" \
    --genome "$genome" \
    --peaks "$peaks" \
    --regions-out "$peaks" \
    --read_shift 4 -5 \
    --extend 100 \
    --prefix "$sample" \
    --outdir "$output/01_atacorrect" \
    --cores "$cores"

# 2. Calculate footprint scores from the corrected cleavage signal.
stage="ScoreBigwig"
TOBIAS ScoreBigwig \
    --signal "$output/01_atacorrect/${sample}_corrected.bw" \
    --regions "$peaks" \
    --output "$output/02_footprint/${sample}_footprints.bw" \
    --score footprint \
    --fp-min 20 \
    --fp-max 50 \
    --flank-min 10 \
    --flank-max 30 \
    --extend 100 \
    --cores "$cores"

# 3. Predict binding sites using the bovine genome, JASPAR motifs and footprint scores.
# Motif/binding p-value cutoffs are software thresholds, not study-level FDR estimates.
stage="BINDetect"
TOBIAS BINDetect \
    --signals "$output/02_footprint/${sample}_footprints.bw" \
    --cond-names "$sample" \
    --motifs "$motifs" \
    --genome "$genome" \
    --peaks "$peaks" \
    --motif-pvalue 1e-4 \
    --bound-pvalue 0.001 \
    --cluster-threshold 0.5 \
    --naming name_id \
    --skip-excel \
    --outdir "$output/03_bindetect" \
    --cores "$cores"

set +x
stage="output validation"
test -s "$output/03_bindetect/bindetect_results.txt"
printf 'Completed TOBIAS analysis for %s\n' "$sample" > "$output/DONE"
echo "Completed: $output"
