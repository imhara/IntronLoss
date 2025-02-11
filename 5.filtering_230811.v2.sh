#!/bin/bash

#updated: 2023.08.11 (Hara Jun) 
# Usage: bash add.sample_name.sh [gene.sam_path] [output_path]

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: bash add.sample_name8.sh [gene.sam_path] [output_path]"
    exit 1
fi

gene_sam_path=$1
output_path=$2

# Make output_dir if it doesn't exist
if [ ! -d ${output_path} ]; then
    mkdir ${output_path}
fi

# Function to process each sample
process_sample() {
    file="$1"
    sample=$(basename "$file" | awk -F'.' '{print $1}')

    if [ -z "$sample" ]; then
        echo "Sample name not found for file: $file"
        return
    fi

    temp_file=$(mktemp)

    echo -e "sample\tQNAME\tFLAG\tRNAME\tPOS\tMAOQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ" > "$temp_file"
    awk -v sample="$sample" -v OFS='\t' 'FNR > 1 && $6 != "*" { print sample, $0 }' "$file" >> "$temp_file"
    awk -F'\t' 'NR > 1 { count[$2]++; if (count[$2] == 2) print $2 } END { close(FILENAME); }' "$temp_file" > "${output_path}/${sample}.qnames.txt"
    awk -F'\t' 'NR==FNR { qnames[$1]; next } $2 in qnames' "${output_path}/${sample}.qnames.txt" "$temp_file" > "${output_path}/${sample}.txt"
    awk -F'\t' '{ pos = $9 - $5; if (pos >= -326 && pos <= 326) print }' "${output_path}/${sample}.txt" > "${output_path}/${sample}.filtered.txt"
    awk -F'\t' -v OFS='\t' 'NR > 1 { print $2"/1" }' "${output_path}/${sample}.filtered.txt" > "${output_path}/${sample}.read1.txt"
    awk -F'\t' -v OFS='\t' 'NR > 1 { print $2"/2" }' "${output_path}/${sample}.filtered.txt" > "${output_path}/${sample}.read2.txt"

    # Remove unnecessary files
    rm -f "${output_path}/${sample}.qnames.txt"
    rm -f "${output_path}/${sample}.txt"
    rm -f "${output_path}/${sample}.filtered.txt"
}

# Loop through the gene.sam files
for file in "${gene_sam_path}"/*.gene.tag.sam; do
    process_sample "$file"
done

echo "Processing complete"

