#!/bin/bash

if [ $# -lt 4 ]
then
        echo usage : $0 [original_sam_file] [read_list.txt] [output_dir] [output_prefix]
        echo "by sample"
        exit 1
fi

original_sam_file=$1
read_list_txt=$2
output_dir=$3
output_prefix=$4

if [ ! -d ${output_dir} ]
then
mkdir ${output_dir}
fi

sample=$(basename ${original_sam_file} .sam)

# Convert BAM to SAM only once
samtools view -@ 2000 -h -o ${output_dir}/${sample}.sam ${original_sam_file}

# Create an empty final SAM file with the desired prefix
> ${output_dir}/${sample}.${output_prefix}.sam

# Filter reads and append to the final SAM file
awk 'NR==FNR{reads[$1]; next} $1 in reads' ${read_list_txt} ${output_dir}/${sample}.sam >> ${output_dir}/${sample}.${output_prefix}.sam

if [ $? -ne 0 ]
then
    echo "An error occurred when extracting. Please check your input files and try again."
    exit 1
fi

# Remove the intermediate SAM file
rm -rf ${output_dir}/${sample}.sam

echo "Extraction complete! Results saved in ${output_dir}/${sample}.${output_prefix}.sam" 

