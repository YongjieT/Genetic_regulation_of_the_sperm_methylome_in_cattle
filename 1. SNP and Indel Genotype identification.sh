#!/usr/bin/env bash
# Staged fastp/GATK 3.8 commands; run per-sample and cohort steps separately.
set -euo pipefail
: "${SAMPLE:?Set SAMPLE}"

# 1. Quality control with fastp

fastp -i ${SAMPLE}/${SAMPLE}_1.raw.fq.gz \
      -o ${SAMPLE}/${SAMPLE}_1.clean.fq.gz \
      -I ${SAMPLE}/${SAMPLE}_2.raw.fq.gz \
      -O ${SAMPLE}/${SAMPLE}_2.clean.fq.gz \
      --html ${SAMPLE}/${SAMPLE}.html \
      --json ${SAMPLE}/${SAMPLE}.json \
      --thread 16

# 2. Base quality score recalibration
java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
    -nct 8 -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
    -I ${SAMPLE}_dedup.bam \
    -knownSites ARS1.2PlusY_BQSR_v3.vcf \
    --bqsrBAQGapOpenPenalty 45 \
    -o ${SAMPLE}.recal.table

# 3. Per-sample gVCF calling (GATK 3.8).
# Input BAM must already have BQSR applied; alignment/BQSR application is upstream.
java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 \
  -jar GenomeAnalysisTK.jar \
  -T HaplotypeCaller -nct 32 \
  -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
  -I "${SAMPLE}_dedup_recal.bam" \
  -o "${SAMPLE}.g.vcf.gz" -ERC GVCF \
  -variant_index_type LINEAR -variant_index_parameter 128000

# 4. Combine GVCFs (this should run once, not per sample)
java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar GenomeAnalysisTK.jar -T CombineGVCFs \
    -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
    -V gvcf.list \
    -o 125Combined.g.vcf.gz

# 5. Chromosome 1 genotyping example; prepare chromosome-specific combined gVCFs.
java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
    -V Combine_chr1.g.vcf.gz \
    -o Genotype_chr1.vcf.gz

# 6. Extract SNPs and indels; indel-specific filters are not provided.
for CHR in {1..29} X; do
    # Select SNPs
    java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar GenomeAnalysisTK.jar -T SelectVariants \
        -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
        --selectTypeToInclude SNP \
        -V Genotype_chr${CHR}.vcf.gz \
        -o rawSNP_chr${CHR}.vcf.gz

    # Select InDels
    java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar GenomeAnalysisTK.jar -T SelectVariants \
        -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
        --selectTypeToInclude INDEL \
        -V Genotype_chr${CHR}.vcf.gz \
        -o rawInDel_chr${CHR}.vcf.gz
done

# 7. Variant filtration (chromosome-specific)
for CHR in {1..29} X; do
    java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar GenomeAnalysisTK.jar -T VariantFiltration \
        -V rawSNP_chr${CHR}.vcf.gz \
        -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
        -filter "QD < 2.0" -filterName "QD2" \
        -filter "FS > 60.0" -filterName "FS60" \
        -filter "MQ < 40.0" -filterName "MQ40" \
        -filter "SOR > 3.0" -filterName "SOR3" \
        -filter "MQRankSum < -12.5" -filterName "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" -filterName "ReadPosRankSum-8" \
        -o refilter_chr${CHR}_snp.vcf.gz
done

# 8. Select passed variants (chromosome-specific)
for CHR in {1..29} X; do
    java -Xmx80G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar GenomeAnalysisTK.jar -T SelectVariants \
        -R Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
        --excludeFiltered \
        -V refilter_chr${CHR}_snp.vcf.gz \
        -o pass_chr${CHR}_snp.vcf.gz
done
