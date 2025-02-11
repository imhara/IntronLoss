#! /bin/bash -e 

##pipeline: Extraction of reads having intron loss patterns 
##authors: Ha Ra Jun
##Date: 2022.09.20
##Updated: 2024.01.22

if [ $# -lt 4 ]
then
	echo usage: $0 [fastq_dir] [fastq_txt] [read_dir] [output_dir]
	exit 1
fi

fastq_dir=$1
fastq_txt=$2
read_dir=$3
output_dir=$4

cat $fastq_txt | while read line;
do 
sample=$line
fastq1=$(echo ${sample} | awk '{print $1}')
fastq2=$(echo ${sample} | awk '{print $2}')
output_prefix=""

if [ ${#sample} -eq 8 ]; then
	output_prefix=${sample%%.R3*} 
else
	output_prefix=${sample%%.R3*} 
fi

echo "$output_prefix" 

#Program
seqkit=/data/data1/Hara/seqkit

#make output_dir

if [ ! -d ${output_dir} ]
then
mkdir ${output_dir}
fi

#Extraction of reads1 

$seqkit grep -f ${read_dir}/${output_prefix}.read1.txt ${fastq_dir}/${fastq1} -o ${output_dir}/${output_prefix}.filtered.R3.r1.fastq 

#Extraction of read2 

$seqkit grep -f ${read_dir}/${output_prefix}.read2.txt ${fastq_dir}/${fastq2} -o ${output_dir}/${output_prefix}.filtered.R3.r2.fastq 



done 
