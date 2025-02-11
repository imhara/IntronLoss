#!/bin/bash -e

##pipeline: fastq_to_maf_CollectHsMetric_GATK4_mutect2_AF0_GGS_true_filtering-new.sh
##authors:Ji Hwe Oh, Chang Ohk Sung, Keyoung-Woon Choi, and Sung-Min Chun
## ALL-IN-ONE pipline GATK4 Mutect2           						-- updated at 2020.10.21
## ADD MFRL info / fix whiteblackanno /filtering2-2(not weight) / MFRL comma to slash 	-- updated at 2020.09.24
## ADD VCF Sorting, Decomposing and Normalizing for multi allele variant              	-- updated at 2020.11.02
## sample names setting 								-- updataed at 2024.01.22



if [ $# -lt 8 ]
then
        echo usage: $0 [fastq_dir] [fastq_txt] [thread] [output_dir] [final_bam_dir] [ref_dir] [exon_ref] [gene_ref]
        exit 1
fi

fastq_dir=$1
fastq_txt=$2
thread=$3
output_dir=$4
final_bam_dir=$5
ref_dir=$6
exon_ref=$7
gene_ref=$8

#export GATK_LOCAL_JAR=/home/tez/anaconda3/pkgs/gatk4-4.1.3.0-0/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar

cat $fastq_txt|while read line;
do
sample=$line
fastq1=$(echo ${sample}|awk '{print $1}')
fastq2=$(echo ${sample}|awk '{print $2}')
output_prefix=""
panel=$(echo ${exon_ref} | cut -d '.' -f 1)
if [ ${#sample} -eq 8 ]; then
	output_prefix=${sample%%.filtered*} 
else
	output_prefix=${sample%%.filtered*}
fi


##programs (anaconda3)

# bwa=/home/tez/anaconda3/bin                                 #0.7.17
# samtools=/home/tez/anaconda3/bin                            #1.9
# picard=/home/tez/anaconda3/bin                              #2.20.5
# gatk4=/home/tez/anaconda3/bin                               #4.1.3.0
gatk4=gatk                                                    
#Indelocator=/home/tez/Programs/GenomeAnalysisTK-2.3-9        
#mutect=/home/tez/Programs/muTect-1.A1.7-bin
#JAVA=/home/tez/anaconda3/bin                                 #8.45.14
#vcf2maf=/home/tez/anaconda3/bin                              #1.6.17


##reference data

ref=~/reference/b37/human_g1k_v37.fasta
mills=~/reference/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
phase1indel=~/reference/b37/1000G_phase1.indels.b37.vcf
dbsnp=~/reference/b37/dbsnp_138.b37.vcf
cosmic=~/reference/b37/Cosmic.hg19.vcf
PON=~/PON/OPv2_mutect_PON.vcf
gnomad=~/reference/b37/af-only-gnomad.raw.sites.b37.vcf


## make output_dir directory

if [ ! -d ${output_dir} ]
then
mkdir ${output_dir}
fi


if [ ! -d ${output_dir}/${output_prefix} ]
then
mkdir ${output_dir}/${output_prefix}
fi

if [ ! -d ${final_bam_dir} ]
then
mkdir ${final_bam_dir}
fi

if [ ! -d ${final_bam_dir}/${panel}.il.filtered.final.sam ]
then
mkdir ${final_bam_dir}/${panel}.il.filtered.final.sam
fi

#if [ ! -d ${output_dir}/${panel}.il.filtered.final.exon.sam ]
#then
#mkdir ${output_dir}/${panel}.il.filtered.final.exon.sam 
#fi

#if [ ! -d ${output_dir}/${panel}.il.filtered.final.gene.sam ]
#then
#mkdir ${output_dir}/${panel}.il.filtered.final.gene.sam
#fi

#export PATH=$HOME/Programs/ensemble-vep-release-88:$PATH


echo "start ${output_prefix} running"
date | tee -a $output_dir/$output_prefix/${output_prefix}_log.txt

RG='@RG\tID:'$output_prefix'\tPL:illumina\tSM:'$output_prefix'\tLB:'$output_prefix''

bwa mem -M -R $RG $ref $fastq_dir/$fastq1 $fastq_dir/$fastq2 -t $thread > $output_dir/$output_prefix/${output_prefix}.${panel}.R3.sam | tee -a $output_dir/$output_prefix/${output_prefix}_log.txt
date | tee -a $output_dir/$output_prefix/${output_prefix}.${panel}_log.txt

echo "convert sam to bam"

samtools view -bhS $output_dir/$output_prefix/${output_prefix}.${panel}.R3.sam -o $output_dir/$output_prefix/${output_prefix}.${panel}.R3.bam | tee -a $output_dir/$output_prefix/${output_prefix}.${panel}_log.txt
date | tee -a $output_dir/$output_prefix/${output_prefix}.${panel}_log.txt


## GATK4 sort >> recalibration 1 2
date
echo "sort reads by coordinate"

$gatk4 --java-options "-Xmx8g"  SortSam \
        --INPUT ${output_dir}/${output_prefix}/${output_prefix}.${panel}.R3.bam \
        --OUTPUT ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.R3.bam \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY SILENT 

echo "remove duplicates"
#$gatk4 --java-options "-Xmx8g"  MarkDuplicates \
#        --INPUT ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.R3.bam \
#        --METRICS_FILE ${output_dir}/${output_prefix}/${output_prefix}.${panel}_R.duplicates \
#        --OUTPUT ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.R3.bam \
#        --VALIDATION_STRINGENCY SILENT 

date | tee -a $output_dir/$output_prefix/${output_prefix}.${panel}_log.txt

echo "recalibration step 1"
#$gatk4 --java-options "-Xmx8g" BaseRecalibrator \
#        --reference ${ref} \
#        --input ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.R3.bam \
#        --known-sites ${dbsnp} \
#        --known-sites ${mills} \
#        --known-sites ${phase1indel} \
#        --output ${output_dir}/${output_prefix}/${output_prefix}.${panel}.R3.recal.table 

date | tee -a $output_dir/$output_prefix/${output_prefix}.${panel}_log.txt

echo "recalibration step 2"
#$gatk4 --java-options "-Xmx8g" ApplyBQSR \
#        --reference ${ref} \
#        --input ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.R3.bam \
#        --bqsr-recal-file ${output_dir}/${output_prefix}/${output_prefix}.${panel}.R3.recal.table \
#        --output ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bam 



samtools index ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.R3.bam ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.R3.bam.bai
#samtools index ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.R3.bam ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.R3.bam.bai 
#samtools index ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bam ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bam.bai


echo "${prefix} bam move"

mv ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.R3.bam ${final_bam_dir}
mv ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.R3.bam.bai ${final_bam_dir}


#mv ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bam ${final_bam_dir}
#mv ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bam.bai ${final_bam_dir}

date | tee -a $output_dir/$output_prefix/${output_prefix}.${panel}_log.txt

#final.sam#
samtools view ${final_bam_dir}/${output_prefix}.${panel}.sorted.R3.bam > ${final_bam_dir}/${output_prefix}.${panel}.sorted.R3.sam
grep -v '^@' ${final_bam_dir}/${output_prefix}.${panel}.sorted.R3.sam > ${final_bam_dir}/${output_prefix}.${panel}.tmp.R3.sam


#samtools view ${final_bam_dir}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bam > ${final_bam_dir}/${output_prefix}.${panel}.sorted.dedup.recal.R3.sam
#grep -v '^@' ${final_bam_dir}/${output_prefix}.${panel}.sorted.dedup.recal.R3.sam > ${final_bam_dir}/${output_prefix}.${panel}.tmp.R3.sam 


awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ${final_bam_dir}/${output_prefix}.${panel}.tmp.R3.sam > ${final_bam_dir}/${output_prefix}.${panel}.il.filtered.final.R3.sam 
mv ${final_bam_dir}/${output_prefix}.${panel}.il.filtered.final.R3.sam ${final_bam_dir}/${panel}.il.filtered.final.sam 
awk '{print "'"$output_prefix"'""\t"$0}' ${final_bam_dir}/${panel}.il.filtered.final.sam/${output_prefix}.${panel}.il.filtered.final.R3.sam > ${final_bam_dir}/${panel}.il.filtered.final.sam/${output_prefix}.${panel}.il.filtered.final.R3.sam2
awk '$10 > 0 {print $0}' ${final_bam_dir}/${panel}.il.filtered.final.sam/${output_prefix}.${panel}.il.filtered.final.R3.sam2 > ${final_bam_dir}/${panel}.il.filtered.final.sam/${output_prefix}.${panel}.R3.r1.sam 
awk '$10 <= 0 {print $0}' ${final_bam_dir}/${panel}.il.filtered.final.sam/${output_prefix}.${panel}.il.filtered.final.R3.sam2 > ${final_bam_dir}/${panel}.il.filtered.final.sam/${output_prefix}.${panel}.R3.r2.sam 

#rm -rf ${output_dir}/${output_prefix}/${output_prefix}.${panel}.sorted.dedup.recal.R3.bai 
#rm -rf ${output_dir}/${output_prefix}/${output_prefix}.${panel}.R3.sam 
##m -rf ${output_dir}/${output_prefix}/*log.txt 
#rm -rf ${final_bam_dir}/${output_prefix}.${panel}.tmp.R3.sam 
#rm -rf ${output_dir}
done
