#!/usr/bin/env bash
# Multi-caller SV discovery and cohort-wide BayesTyper genotyping.
# Usage: bash SV_calling_joint_BayesTyper.sh call|prepare|candidate|kmc|genotype|merge|qc [batch_samples.txt]
# Run stages separately; kmc and genotype can use multiple sample batches.
set -euo pipefail

stage=${1:-}
batch_file=${2:-}
case "$stage" in
    call|prepare|candidate|kmc|genotype|merge|qc) ;;
    *) echo "Usage: $0 call|prepare|candidate|kmc|genotype|merge|qc [batch_samples.txt]" >&2; exit 1 ;;
esac

# Configure inputs and working directories before running.
x=$(cat samplelist.txt)
reference=/path/to/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
bam_dir=/path/to/ALLbam
work=/path/to/SV_work
upstream=/path/to/prepared_1MBlength_vcfs
filtered=/path/to/existing_norm_filtered_vcfs
ploidy=/path/to/cattle_ploidy.tsv
kmc_bin=kmc
threads=${SLURM_CPUS_PER_TASK:-6}

source /path/to/miniconda3/etc/profile.d/conda.sh
mkdir -p "$work"

# 1. Identify SVs independently in each sample with five callers.
if [[ "$stage" == call ]]; then
    for Bam in $x
    do
        bam="$bam_dir/newRG_${Bam}_deup.bam"
        dir="$work/Calling/$Bam"
        [[ ! -e "$dir" ]] || { echo "Output exists: $dir" >&2; exit 1; }
        mkdir -p "$dir/Dysgu_tmp"

        conda activate SV
        configManta.py --bam "$bam" --referenceFasta "$reference" --runDir "$dir/manta"
        python "$dir/manta/runWorkflow.py" -j "$threads"
        bcftools view -Ov -o "$dir/${Bam}_manta.vcf" \
            "$dir/manta/results/variants/diploidSV.vcf.gz"

        delly call -g "$reference" -o "$dir/${Bam}_delly.bcf" "$bam"
        bcftools view -Ov -o "$dir/${Bam}_delly.vcf" "$dir/${Bam}_delly.bcf"

        whamg -x 6 -a "$reference" -f "$bam" \
            > "$dir/${Bam}_Wham.vcf" 2> "$dir/Wham.err"

        smoove call --outdir "$dir/smoove" --name "$Bam" \
            --fasta "$reference" -p "$threads" --genotype "$bam"
        bcftools view -Ov -o "$dir/${Bam}_smoove.vcf" \
            "$dir/smoove/${Bam}-smoove.genotyped.vcf.gz"

        conda activate Dysgu
        dysgu run -p "$threads" "$reference" "$dir/Dysgu_tmp" "$bam" \
            > "$dir/${Bam}_Dysgu.vcf"
    done
fi

# 2. Convert and normalize alleles using the existing upstream size-filtered VCFs.
# Existing length/quality filtering rules are unchanged and are NOT redefined here.
if [[ "$stage" == prepare ]]; then
    conda activate WGS3
    for SampleID in $x
    do
        dir="$work/Prepared/$SampleID"
        [[ ! -e "$dir" ]] || { echo "Output exists: $dir" >&2; exit 1; }
        mkdir -p "$dir"
        for caller in delly manta smoove Wham Dysgu
        do
            input="$upstream/$SampleID/${SampleID}_${caller}_1MBlength.vcf"
            prefix="$dir/${SampleID}_${caller}_convert"
            bayesTyperTools convertAllele -v "$input" -g "$reference" -z -o "$prefix"
            bcftools norm -c x -f "$reference" -m -any -w 1000000 \
                -Ov -o "$dir/${SampleID}_${caller}_norm.vcf" "$prefix.vcf.gz"
        done
    done
fi

# 3. Combine SV candidates from ALL samples and ALL five callers.
# Input *_norm_filtered.vcf files are generated using the unchanged prior filters.
if [[ "$stage" == candidate ]]; then
    conda activate WGS3
    dir="$work/Candidate"
    [[ ! -e "$dir" ]] || { echo "Output exists: $dir" >&2; exit 1; }
    mkdir -p "$dir/input"
    chromosomes="$(seq -s, 1 29),X"
    inputs=()
    for SampleID in $x
    do
        for entry in Delly:delly Dysgu:DysguSorted Manta:manta Smoove:smoove Wham:Wham
        do
            caller=${entry%%:*}
            suffix=${entry#*:}
            vcf="$filtered/$SampleID/${SampleID}_${suffix}_norm_filtered.vcf"
            [[ -s "$vcf" ]] || { echo "Missing VCF: $vcf" >&2; exit 1; }
            selected="$dir/input/${SampleID}_${caller}.vcf"
            bcftools view -t "$chromosomes" -Ov -o "$selected" "$vcf"
            inputs+=("$caller:$selected")
        done
    done
    combined=$(IFS=,; echo "${inputs[*]}")
    bayesTyperTools combine -v "$combined" -o "$dir/cohort" -z
fi

# 4. Count 55-mers and generate Bloom filters for each sample.
# Optional batch_samples.txt lists sample IDs, one per line; batches are independent.
if [[ "$stage" == kmc ]]; then
    conda activate WGS3
    ids=$x
    [[ -z "$batch_file" ]] || ids=$(cat "$batch_file")
    # KMC BAM input must contain unmapped reads as well as mapped reads.
    for SampleID in $ids
    do
        dir="$work/KMC/$SampleID"
        prefix="$dir/$SampleID"
        [[ ! -e "$dir" ]] || { echo "Output exists: $dir" >&2; exit 1; }
        mkdir -p "$dir"
        "$kmc_bin" -k55 -t"$threads" -m50 -ci1 -fbam \
            "$bam_dir/${SampleID}_deup.bam" "$prefix" "$dir"
        bayesTyperTools makeBloom -k "$prefix" -p "$threads"
    done
fi

# 5. Genotype each batch against the SAME cohort-wide SV candidate VCF.
# batch_samples.txt contains sample IDs only; all bulls are male.
if [[ "$stage" == genotype ]]; then
    conda activate WGS3
    [[ -s "$batch_file" ]] || { echo 'Provide batch_samples.txt' >&2; exit 1; }
    batch=$(basename "$batch_file" .txt)
    dir="$work/Genotype/$batch"
    [[ ! -e "$dir" ]] || { echo "Output exists: $dir" >&2; exit 1; }
    mkdir -p "$dir"
    : > "$dir/samples.tsv"
    for SampleID in $(cat "$batch_file")
    do
        grep -Fxq -- "$SampleID" samplelist.txt || { echo "Unknown sample: $SampleID" >&2; exit 1; }
        prefix="$work/KMC/$SampleID/$SampleID"
        for ext in kmc_pre kmc_suf bloomMeta bloomData; do
            [[ -s "$prefix.$ext" ]] || { echo "Missing KMC input: $prefix.$ext" >&2; exit 1; }
        done
        printf '%s\tM\t%s\n' "$SampleID" "$prefix" >> "$dir/samples.tsv"
    done
    candidate="$work/Candidate/cohort.vcf.gz"
    [[ -s "$candidate" && -s "$ploidy" ]] || { echo 'Missing cohort VCF or ploidy file' >&2; exit 1; }
    cd "$dir"
    bayesTyper cluster -v "$candidate" -s samples.tsv -g "$reference" -p "$threads"
    mapfile -t units < <(find "$dir" -maxdepth 1 -type d -name 'bayestyper_unit_*' | sort -V)
    [[ ${#units[@]} -gt 0 ]] || { echo 'No BayesTyper units' >&2; exit 1; }
    vcfs=()
    for unit in "${units[@]}"
    do
        bayesTyper genotype -v "$unit/variant_clusters.bin" \
            -c "$dir/bayestyper_cluster_data" -s "$dir/samples.tsv" \
            -g "$reference" --chromosome-ploidy-file "$ploidy" \
            -o "$unit/bayestyper" -z -p "$threads"
        vcfs+=("$unit/bayestyper.vcf.gz")
    done
    bcftools concat -Oz -o "$dir/${batch}.vcf.gz" "${vcfs[@]}"
    tabix -p vcf "$dir/${batch}.vcf.gz"
fi

# 6. Merge genotyped batches into one population VCF.
# Four-column manifest: sample ID, source, batch, absolute VCF path.
# Choose one source VCF per sample if a sample occurs in multiple batches.
if [[ "$stage" == merge ]]; then
    conda activate WGS3
    manifest="$work/selected_sample_sources.tsv"
    [[ -s "$manifest" ]] || { echo "Missing merge manifest: $manifest" >&2; exit 1; }
    dir="$work/Merge"
    [[ ! -e "$dir" ]] || { echo "Output exists: $dir" >&2; exit 1; }
    diff <(sort -u samplelist.txt) <(cut -f1 "$manifest" | sort -u)
    [[ $(cut -f1 "$manifest" | wc -l) -eq $(cut -f1 "$manifest" | sort -u | wc -l) ]] || { echo "Duplicate samples in manifest" >&2; exit 1; }
    mkdir -p "$dir"
    vcfs=()
    index=0
    signature=''
    while IFS= read -r input
    do
        index=$((index + 1))
        keep="$dir/keep_$index.txt"
        subset="$dir/input_$index.vcf.gz"
        awk -F '\t' -v p="$input" '$4==p {print $1}' "$manifest" > "$keep"
        bcftools view -S "$keep" -Oz -o "$subset" "$input"
        tabix -p vcf "$subset"
        current=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$subset" | md5sum | cut -d' ' -f1)
        [[ -z "$signature" || "$signature" == "$current" ]] || { echo "Batch variant sets differ" >&2; exit 1; }
        signature=$current
        vcfs+=("$subset")
    done < <(cut -f4 "$manifest" | sort -u)
    bcftools merge --threads "$threads" --filter-logic x --info-rules ACP:max \
        -Oz -o "$dir/RAW.vcf.gz" "${vcfs[@]}"
    bcftools view -f PASS -Oz -o "$dir/PASS.vcf.gz" "$dir/RAW.vcf.gz"
    tabix -p vcf "$dir/PASS.vcf.gz"
    diff <(sort samplelist.txt) <(bcftools query -l "$dir/PASS.vcf.gz" | sort)
fi
