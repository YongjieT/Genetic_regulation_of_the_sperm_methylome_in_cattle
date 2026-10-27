#!/usr/bin/env bash
# Staged STR commands; configure each working directory and caller path.
set -euo pipefail

# 1. GangSTR Analysis

# 1.1 Run GangSTR
source activate STR
cd gangSTR

OUT=gangSTR
Ref=Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
Ref_trf=Bos_taurus.GangSTR.bed

i=${SAMPLE:?Set SAMPLE}

GangSTR --bam newRG_${SAMPLE}_deup.bam \
        --ref ${Ref} \
        --regions ${Ref_trf} \
        --out ${OUT}/${i}_GangSTR_calls

# 1.2 Process VCF files
cd gangSTR
for vcf in *.vcf; do
    [[ -f "$vcf" ]] || { echo "ERROR: No VCF inputs" >&2; exit 1; }
    bcftools sort "$vcf" -Oz -o "${vcf}.gz"
    tabix -p vcf "${vcf}.gz"
done
ls *.vcf.gz > 125Bull_STR_list.txt

# 1.3 Merge VCF files
cd gangSTR/2Merge
conda activate STR
mergeSTR \
  --vcfs-list 125Bull_STR_list.txt \
  --out 125bull_STR_gangstr \
  --vcftype gangstr

zcat 125bull_STR_gangstr.vcf.gz | vcf-sort | bgzip -c > 125bull_STR_gangstr.sorted.vcf.gz
tabix -p vcf 125bull_STR_gangstr.sorted.vcf.gz

# 1.4 Filter STR calls
dumpSTR \
    --gangstr-filter-spanbound-only \
    --gangstr-filter-badCI \
    --gangstr-max-call-DP 1000 \
    --gangstr-min-call-DP 10 \
    --min-locus-callrate 0.8 \
    --min-locus-hwep 0.000001 \
    --drop-filtered \
    --vcftype gangstr \
    --zip \
    --vcf 125bull_STR_gangstr.sorted.vcf.gz \
    --out 125bull_STR_gangstr_filter

# 2. ExpansionHunter Analysis

ExpansionHunter=ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter
cd ExpansionHunter

OUT=ExpansionHunter
Ref=Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
Ref_tf=Reference/Bos_taurus.STR_fot_EH.json

# Run ExpansionHunter
    $ExpansionHunter --reads newRG_${SAMPLE}_deup.bam \
                    --reference ${Ref} \
                    --threads 4 \
                    --variant-catalog ${Ref_tf} \
                    --output-prefix ExpansionHunter/${SAMPLE}_ExpansionHunter_calls

# Process VCF files
cd ExpansionHunter
for vcf in *.vcf; do
    [[ -f "$vcf" ]] || { echo "ERROR: No VCF inputs" >&2; exit 1; }
    bcftools sort "$vcf" -Oz -o "${vcf}.gz"
    tabix -p vcf "${vcf}.gz"
done
ls *.vcf.gz > 125Bull_STR_list.txt

# Merge and filter ExpansionHunter results
mergeSTR \
  --vcfs-list 125Bull_STR_list.txt \
  --out 125bull_STR_EH \
  --vcftype eh

zcat 125bull_STR_EH.vcf.gz | vcf-sort | bgzip -c > 125bull_STR_EH.sorted.vcf.gz
tabix -p vcf 125bull_STR_EH.sorted.vcf.gz

dumpSTR \
    --eh-min-call-LC 10 \
    --min-locus-callrate 0.8 \
    --min-locus-hwep 0.000001 \
    --vcftype eh \
    --drop-filtered \
    --zip \
    --vcf 125bull_STR_EH.sorted.vcf.gz \
    --out 125bull_STR_EH_filter

# 3. HipSTR Analysis

cd HipSTR
OUT=HipSTR
Ref=Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
Ref_trf=Reference/Bos_taurus.HipSTR.bed

# Run HipSTR
HipSTR=HipSTR/HipSTR
$HipSTR --bam-files SampleList125Bull.txt \
        --fasta ${Ref} \
        --regions ${Ref_trf} \
        --str-vcf ${OUT}/125Bull_HipSTR_calls.vcf.gz \
        --max-str-len 200

# Filter HipSTR results
dumpSTR \
    --hipstr-max-call-DP 1000 \
    --hipstr-min-call-DP 10 \
    --min-locus-callrate 0.8 \
    --min-locus-hwep 0.000001 \
    --vcftype hipstr \
    --drop-filtered \
    --zip \
    --vcf 125Bull_HipSTR_calls.vcf.gz \
    --out 125Bull_STR_hipstr_filter

# Apply HipSTR correction
python3 HipSTR/HipSTR_correction/Hipstr_correction.py \
HipSTR/125Bull_STR_hipstr_filter.vcf.gz 125Bull_STR_hipstr_filter_corrected.vcf
bgzip -c 125Bull_STR_hipstr_filter_corrected.vcf > 125Bull_STR_hipstr_filter_corrected.vcf.gz
tabix -p vcf 125Bull_STR_hipstr_filter_corrected.vcf.gz

# 4. EnsembleTR Integration

cd EnsembleTR
source activate EnsembleTR

Results_gangSTR=125bull_STR_gangstr_filter.vcf.gz
Results_HipSTR=125Bull_STR_hipstr_filter.vcf.gz
Results_EH=125bull_STR_EH.sorted_filter.modified.vcf.gz

# Supply the modified EH VCF and choose the intended corrected/uncorrected HipSTR input.
# Run EnsembleTR
ensembletr --out Hip_Gang_EH.vcf \
--ref Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
--vcfs ${Results_gangSTR},${Results_HipSTR},${Results_EH}

# Confirm that the integrated VCF is compatible with --vcftype hipstr.
# Final filtering
dumpSTR \
    --min-locus-callrate 0.8 \
    --vcftype hipstr \
    --zip \
    --vcf EnsembleTR/Hip_Gang_EH.vcf \
    --out Hip_Gang_EH_callrate0.8
