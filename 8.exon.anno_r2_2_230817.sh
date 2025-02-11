#! /bin/bash -e

#pipeline
#Hara Jun
#2023.03.28
#Updated: 2024.01.22

if [ $# -lt 3 ]
then
        echo usage : $0 [gene_sam_dir] [gene_sam_txt] [output_dir]
        exit 1
fi
gene_sam_dir=$1
gene_sam_txt=$2
output_dir=$3

if [ ! -d ${output_dir} ]
then
mkdir ${output_dir} 
fi 


while read data;
do
    echo "data:$data"
    while read line;
    do
        sample_txt=$(awk '{print $1}' <<< $line)
        echo $sample_txt
        sample="" 
	if [ ${#sample_txt} -eq 8 ]; then
		sample=${sample_txt%%.op*}
	else
		sample=${sample_txt%%.op*} 
	fi 
	echo $sample
	read=$(awk '{print $2}' <<< $line) 
	echo $read 
        gene=$(awk '{print $4}' <<< $line)
        echo $gene
        pos=$(awk '{print $5}' <<< $line)
	pos2=$(awk '{print $5}' <<< $line) 
        echo $pos2
	flag=$(awk '{print $3}' <<< $line) 
	maoq=$(awk '{print $6}' <<< $line) 
	cigar=$(awk '{print $7}' <<< $line) 
	rnext=$(awk '{print $8}' <<< $line) 
	pnext=$(awk '{print $9}' <<< $line)
	tlen=$(awk '{print $10}' <<< $line)
	seq=$(awk '{print $11}' <<< $line) 


	
awk -v gene_name="$gene" -v position="$pos2" '$1==gene_name && $2<=position && $3>=position  {print}'  /data/data1/Hara/Intron_Loss/2023/target/op_amc_v3_exon.bed.txt | awk -v sample="$sample" -v read="$read" -v flag="$flag" -v gene="$gene" -v pos="$pos" -v maoq="$maoq" -v cigar="$cigar" -v rnext="$rnext" -v pnext="$pnext" -v tlen="$tlen" -v seq="$seq" '{print sample"\t"read"\t"flag"\t"gene"\t"pos"\t"maoq"\t"cigar"\t"rnext"\t"pnext"\t"tlen"\t"seq"\t""r2""\t"$0}' >> ${output_dir}/${sample}.R3.r2.out2.txt

 


    done < ${gene_sam_dir}/$data
done < $gene_sam_txt
