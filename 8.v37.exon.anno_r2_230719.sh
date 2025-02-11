#! /bin/bash -e

#pipeline
#Hara Jun
#2023.04.06
#2023.07.12 updated in super7 server 
#2023.07.17 updated in super7 server 
## replace exon.bed to op_amc_v3.exon3.bed 


if [ $# -lt 3 ]
then
        echo usage : $0 [b37_sam_dir] [b37_sam_txt] [output_dir]
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
    echo $data
    while read line;
    do
        sample_txt=$(awk '{print $1}' <<< $line)
        #echo $sample_txt
	sample="" 
	if [ ${#sample_txt} -eq 8 ]; then
		sample=${sample_txt%%.op*}
	else
		sample=${sample_txt%%.op*} 
	fi 
        read=$(awk '{print $2}' <<< $line)
	chr=$(awk '{print $4}' <<< $line)
        pos=$(awk '{print $5}' <<< $line)
        flag=$(awk '{print $3}' <<< $line)
        maoq=$(awk '{print $6}' <<< $line)
        cigar=$(awk '{print $7}' <<< $line)
        rnext=$(awk '{print $8}' <<< $line)
        pnext=$(awk '{print $9}' <<< $line)
        tlen=$(awk '{print $10}' <<< $line)
        seq=$(awk '{print $11}' <<< $line)

        bed_info=$(awk -v chromosome="$chr" -v position="$pos" '$1==chromosome && $2<=position && $3>=position  {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7}' /data/data1/Hara/Intron_Loss/2023/target/op_amc_v3_v37.bed.txt)


if [ -z "$bed_info" ]
then
    continue
fi
        echo -e "$sample\t$read\t$flag\t$chr\t$pos\t$maoq\t$cigar\t$rnext\t$pnext\t$tlen\t$seq\t$bed_info\tr2" | tr -d '\r' >> ${output_dir}/${sample}.b37.R3.r2.out.txt

    done < ${gene_sam_dir}/$data
done < $gene_sam_txt

