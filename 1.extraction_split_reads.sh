#!/bin/bash -e

# 입력 파라미터 확인
if [ $# -lt 4 ]; then
    echo "Usage: $0 [fastq_dir] [fastq_txt] [thread] [output_dir]"
    exit 1
fi

fastq_dir=$1
fastq_txt=$2
thread=$3
output_dir=$4

# 도구 및 참조 데이터
bwa="bwa"
seqkit=/data/data1/Hara/seqkit
ref=~/reference/b37/human_g1k_v37.fasta

# 출력 디렉토리 생성
mkdir -p ${output_dir}/Rdat/tmp.sam
mkdir -p ${output_dir}/Rdat/v37.sam
mkdir -p ${output_dir}/R3.fastq
mkdir -p ${output_dir}/log

# 분석 시작 시간
start_time=$(date '+%Y-%m-%d %H:%M:%S')
echo "Analysis started at: $start_time"

# Step 1: Mapping 및 SAM 파일 처리
cat ${fastq_txt} | while read sample; do
    fastq1=$(echo ${sample} | awk '{print $1}')
    fastq2=$(echo ${sample} | awk '{print $2}')
    output_prefix=$(echo ${fastq1} | cut -d '.' -f 1)
    timestamp=$(date '+%Y%m%d%H%M%S')
    log_file="${output_dir}/log/${output_prefix}.log_${timestamp}.txt"
    {
        echo "Processing sample: ${output_prefix}"

        # BWA Mapping
        RG="@RG\tID:${output_prefix}\tPL:illumina\tSM:${output_prefix}\tLB:${output_prefix}"
        $bwa mem -M -R "$RG" $ref $fastq_dir/$fastq1 $fastq_dir/$fastq2 -t $thread \
            > ${output_dir}/${output_prefix}.sam
        echo "Generated SAM file: ${output_dir}/${output_prefix}.sam"

        # SAM column ordering
        if [ -f ${output_dir}/${output_prefix}.sam ]; then
            grep -v '^@' ${output_dir}/${output_prefix}.sam \
                > ${output_dir}/Rdat/tmp.sam/${output_prefix}.tmp.sam
            echo "Generated tmp.sam: ${output_dir}/Rdat/tmp.sam/${output_prefix}.tmp.sam"

            awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' \
                ${output_dir}/Rdat/tmp.sam/${output_prefix}.tmp.sam \
                > ${output_dir}/Rdat/v37.sam/${output_prefix}.v37.R.sam
        else
            echo "Error: SAM file not found for ${output_prefix}"
        fi

        # Step 2: R3 스플릿 리드 필터링
        Rscript /data/data1/Hara/Intron_Loss/2024/script/2023/strategy3_v2/2.filtering_2.R \
            ${output_prefix} ${output_dir}/Rdat/v37.sam ${output_dir}
        echo "R3 filtering completed for: ${output_prefix}"
	
	# step 4: read1.txt & read2.txt file merge
	cat ${output_dir}/*_r1.txt > ${output_dir}/read1.txt
	cat ${output_dir}/*_r2.txt > ${output_dir}/read2.txt

        # Step 5: SeqKit을 사용한 리드 추출
        echo "Extracting reads for sample: ${output_prefix}"
        $seqkit grep -f ${output_dir}/read1.txt ${fastq_dir}/${fastq1} \
            -o ${output_dir}/R3.fastq/${output_prefix}.R3.r1.fastq
        echo "Generated R3.r1.fastq for ${output_prefix}"

        $seqkit grep -f ${output_dir}/read2.txt ${fastq_dir}/${fastq2} \
            -o ${output_dir}/R3.fastq/${output_prefix}.R3.r2.fastq
        echo "Generated R3.r2.fastq for ${output_prefix}"
    } 2>&1 | tee -a $log_file

done

# Step 5: 불필요한 파일 정리
{
    echo "Cleaning up temporary and intermediate files..."

    # *.sam 파일 삭제
    find ${output_dir} -name "*.sam" -type f -delete
    echo "Deleted all SAM files."

    # *.txt 파일 tmp 디렉토리로 이동
    mkdir -p ${output_dir}/tmp
    find ${output_dir} -maxdepth 1 -name "*.txt" -exec mv {} ${output_dir}/tmp/ \;
    echo "Moved all TXT files to tmp directory."
} 2>&1 | tee -a ${output_dir}/log/cleanup.log_$(date '+%Y%m%d%H%M%S').txt

# Step 6: fastq_txt 생성

cd ${output_dir}/R3.fastq
ls *r1* | xargs -n 1 basename > r1.txt
ls *r2* | xargs -n 1 basename > r2.txt
paste r1.txt r2.txt > r3.fastq.list.txt
rm -rf ${output_dir}/R3.fastq/r1.txt ${output_dir}/R3.fastq/r2.txt 

# 분석 마친 시간
end_time=$(date '+%Y-%m-%d %H:%M:%S')
echo "Analysis completed at: $end_time"

