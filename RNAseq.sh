#!/bin/bash
fastqc *.fastq
multiqc .
hisat2-build Bos_taurus.fa Bos_taurus.idx
for i in  `ls *.fastq | sed 's/_[12].fastq//g' |sort -u`
do
	hisat2 -x Bos_taurus.idx -1 ${i}_1.fastq -2 ${i}_2.fastq -S ${i}.sam
	picard SortSam -INPUT ${i}.sam -OUTPUT ${i}_sorted.bam -SORT_ORDER coordinate
	samtools index ${i}_sorted.bam
	picard MarkDuplicates -INPUT ${i}_sorted.bam -OUTPUT ${i}_sorted_dedup.bam -METRICS_FILE dedup_metrics.txt
	picard AddOrReplaceReadGroups -I ${i}_sorted_dedup.bam -O ${i}_sorted_dedup_RG.bam -RGID 1 -RGLB lib2 -RGPL illumina -RGPU unit1 -RGSM 3
	picard BuildBamIndex -INPUT ${i}_sorted_dedup_RG.bam
done

gatk CreateSequenceDictionary -R Bos_taurus.fa
samtools faidx Bos_taurus.fa
wget https://ftp.ensembl.org/pub/release-108/variation/vcf/bos_taurus/bos_taurus.vcf.gz
gzip -d bos_taurus.vcf.gz
bgzip bos_taurus.vcf
tabix -f -p vcf bos_taurus.vcf.gz

for f in `ls *.fastq | sed 's/_[12].fastq//g' |sort -u`
do
	gatk SplitNCigarReads -R Bos_taurus.fa -I ${f}_sorted_dedup_RG.bam -O ${f}_sorted_dedup_RG_splitreads.bam
	gatk BaseRecalibrator -I ${f}_sorted_dedup_RG_splitreads.bam -R Bos_taurus.fa --known-sites bos_taurus.vcf.gz -O recal_data.table
	gatk ApplyBQSR -R Bos_taurus.fa -I ${f}_sorted_dedup_RG_splitreads.bam --bqsr-recal-file recal_data.table -O ${f}_recal_reads.bam
	gatk HaplotypeCaller -R Bos_taurus.fa -I ${f}_recal_reads.bam -O ${f}.g.vcf.gz -ERC GVCF
	#gatk CombineGVCFs -R Bos_taurus.fa --variant ${f}.g.vcf.gz -O ${f}.g.vcf.gz
	gatk GenotypeGVCFs -R Bos_taurus.fa -V ${f}.g.vcf.gz -O ${f}_raw_variants.vcf.gz
	gzip -d ${f}_raw_variants.vcf.gz
	gatk SelectVariants -R Bos_taurus.fa -V ${f}_raw_variants.vcf --select-type-to-include SNP -O ${f}_raw_snps.vcf
	gatk SelectVariants -R Bos_taurus.fa -V ${f}_raw_variants.vcf --select-type-to-include INDEL -O ${f}_raw_indels.vcf
	gatk VariantFiltration -V ${f}_raw_snps.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter"FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum <-12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum <-8.0" --filter-name "ReadPosRankSum-8" -O ${f}_snps_filtered.vcf
	gatk VariantFiltration -V ${f}_raw_indels.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum <-20.0" --filter-name "ReadPosRankSum-20" -O ${f}_indels_filtered.vcf
done
