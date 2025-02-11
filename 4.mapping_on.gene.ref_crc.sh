#! /bin/bash -e

##pipeline: fastq_to_sam with GNA11.target.Exon.region.fasta 
##outhors: Ha Ra Jun 
##Date: 2023.01.19
##Updated Date: 2024.01.22

##TEST to find intron loss patterns 

if [ $# -lt 5 ]
then 
	echo usage: $0 [fastq_dir] [fastq_txt] [output_dir] [ref_dir] [gene_ref]
	exit 1
fi

fastq_dir=$1
fastq_txt=$2
output_dir=$3
ref_dir=$4
gene_ref=$5


export GATK_LOCAL_JAR=/home/tez/anaconda3/pkgs/gatk4-4.1.3.0-0/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar

cat $fastq_txt|while read line;
do
sample=$line
fastq1=$(echo ${sample}| awk '{print $1}')
fastq2=$(echo ${sample}| awk '{print $2}') 
output_prefix=""
panel=$(echo ${gene_ref} | cut -d '.' -f 1) 

if [ ${#sample} -eq 8 ]; then
	output_prefix=${sample%%.r1*}
else
	output_prefix=${sample%%.r1*} 
fi

echo "sample:" $sample

# make output_dir 

if [ ! -d ${output_dir} ] 
then 
mkdir ${output_dir}
fi

if [ ! -d ${output_dir}/${output_prefix} ] 
then
mkdir ${output_dir}/${output_prefix} 
fi

if [ ! -d ${output_dir}/Rdat ]
then 
mkdir ${output_dir}/Rdat 
fi

#if [ ! -d ${output_dir}/Rdat/exon.R.sam ]
#then
#mkdir ${output_dir}/Rdat/exon.R.sam
#fi

if [ ! -d ${output_dir}/Rdat/gene.R.sam ]
then 
mkdir ${output_dir}/Rdat/gene.R.sam 
fi

#if [ ! -d ${output_dir}/Rdat/${panel}.exon.tag.sam ]
#then
#mkdir ${output_dir}/Rdat/${panel}.exon.tag.sam 
#fi

if [ ! -d ${output_dir}/Rdat/${panel}.gene.tag.sam ]
then
mkdir ${output_dir}/Rdat/${panel}.gene.tag.sam 
fi 

#if [ ! -d ${output_dir}/Rdat/exon.R.sam_tmp ]
#then
#mkdir ${output_dir}/Rdat/exon.R.sam_tmp
#fi

if [ ! -d ${output_dir}/Rdat/gene.R.sam_tmp ]
then
mkdir ${output_dir}/Rdat/gene.R.sam_tmp
fi

#if [ ! -d ${output_dir}/Rdat/exon_tmp.sam ]
#then
# mkdir ${output_dir}/Rdat/exon_tmp.sam
#fi

if [ ! -d ${output_dir}/Rdat/gene_tmp.sam ]
then
mkdir ${output_dir}/Rdat/gene_tmp.sam 
fi

echo "start ${output_prefix} Mapping to ${panel}.fasta" 
RG='@RG\tID:'$output_prefix'\tPL:illumina\tSM:'$output_prefix'\tLB:'$output_prefix''

# Mapping with chr.ref.fasta & exon.ref.fasta #
echo ${output_dir} ${output_prfix}.${panel}.exon.sam 
#bwa mem -M I-R $RG $ref_dir/$exon_ref $fastq_dir/$fastq1 $fastq_dir/$fastq2 > ${output_dir}/${output_prefix}/${output_prefix}.${panel}.exon.sam
bwa mem -M -R $RG $ref_dir/$gene_ref $fastq_dir/$fastq1 $fastq_dir/$fastq2 > ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.sam

#Prepatarion for using in R with SAM format#

#grep -v '^@' ${output_dir}/${output_prefix}/${output_prefix}.${panel}.exon.sam > ${output_dir}/Rdat/exon_tmp.sam/${output_prefix}.${panel}.exon_tmp.sam 
grep -v '^@' ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.sam >  ${output_dir}/Rdat/gene_tmp.sam/${output_prefix}.${panel}.gene_tmp.sam

#awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ${output_dir}/Rdat/exon_tmp.sam/${output_prefix}.${panel}.exon_tmp.sam > ${output_dir}/Rdat/exon.R.sam/${output_prefix}.${panel}.exon.R.sam 
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ${output_dir}/Rdat/gene_tmp.sam/${output_prefix}.${panel}.gene_tmp.sam > ${output_dir}/Rdat/gene.R.sam/${output_prefix}.${panel}.gene.R.sam

#SAM to BAM#
echo "Convert .sam to .bam"
#samtools view -bhS ${output_dir}/${output_prefix}/${output_prefix}.${panel}.exon.sam -o ${output_dir}/${output_prefix}/${output_prefix}.${panel}.exon.bam
#samtools view -bhS ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.sam -o ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.bam 

#gatk --java-options "-Xmx8g" SortSam \
#     --INPUT ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.bam \
#     --OUTPUT ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.sorted.bam \
#     --SORT_ORDER coordinate \
#     --VALIDATION_STRINGENCY SILENT
#samtools index ${output_dir}/${output_prefix}/${output_prefix}.${panel}.gene.sorted.bam.bai

#awk '{if($3 != "*" && $7 != "*") print $0}' ${output_dir}/Rdat/exon.R.sam/${output_prefix}.${panel}.exon.R.sam > ${output_dir}/Rdat/exon.R.sam_tmp/${output_prefix}.${panel}.exon.R.sam_tmp 
#awk '{if($7 != "=") print $0}' ${output_dir}/Rdat/exon.R.sam_tmp/${output_prefix}.${panel}.exon.R.sam_tmp > ${output_dir}/Rdat/${panel}.exon.tag.sam/${output_prefix}.${panel}.exon.tag.sam 

awk '{if($3 != "*" && $7 != "*") print $0}' ${output_dir}/Rdat/gene.R.sam/${output_prefix}.${panel}.gene.R.sam > ${output_dir}/Rdat/gene.R.sam_tmp/${output_prefix}.${panel}.gene.R.sam_tmp
awk '{if($7 == "=") print $0}' ${output_dir}/Rdat/gene.R.sam_tmp/${output_prefix}.${panel}.gene.R.sam_tmp > ${output_dir}/Rdat/${panel}.gene.tag.sam/${output_prefix}.${panel}.gene.tag.sam 


#find ${output_dir}/Rdat/${panel}.exon.tag.sam/ -size 0c -delete

rm -rf ${output_dir}/*smp*
#rm -rf ${output_dir}/Rdat/exon.R.sam 
rm -rf ${output_dir}/Rdat/gene.R.sam
#rm -rf ${output_dir}/Rdat/exon.R.sam_tmp
rm -rf ${output_dir}/Rdat/gene.R.sam_tmp
#rm -rf ${output_dir}/Rdat/exon_tmp.sam
rm -rf ${output_dir}/Rdat/gene_tmp.sam 


done
