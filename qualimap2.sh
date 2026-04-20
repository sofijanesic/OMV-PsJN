#!/bin/bash

#SBATCH -J qualimap_igor
#SBATCH -c 30
#SBATCH --mem 48G
#SBATCH -o /lustre/imgge/kurs_april/igor/logs/%j.out
#SBATCH -e /lustre/imgge/kurs_april/igor/logs/%j.err

module load qualimap
#module load gatk
#module load R R-bundle-Bioconductor R-bundle-CRAN

OUTPUT_DIR="/lustre/imgge/kurs_april/igor/3_rna_seq_marija_s/bams/qualimap"
INPUT_DIR="/lustre/imgge/kurs_april/igor/3_rna_seq_marija_s/bams/3_star"
GTF="/lustre/imgge/db/hg38/ensembl_v113/genome.gtf"

let "java_mem = $SLURM_MEM_PER_NODE / 1024 - 4"
export JAVA_OPTS="-Xms4G -Xmx${java_mem}G -XX:MaxMetaspaceSize=4G"

if [ ! -d ${OUTPUT_DIR} ];then
	mkdir ${OUTPUT_DIR}
fi

files=$(find ${INPUT_DIR} -name *Aligned.sortedByCoord.out.bam)


for file in ${files};do
	file_n=$(basename ${file})
	qualimap rnaseq \
		-pe \
		-bam ${INPUT_DIR}/${file_n} \
		-outdir ${OUTPUT_DIR}/${file_n} \
		-gtf ${GTF}
	qualimap bamqc -bam ${INPUT_DIR}/${file_n} -nt 30 -gff $GTF -outdir ${OUTPUT_DIR}/${file_n}_bamqc2
#	gatk CollectInsertSizeMetrics -I ${INPUT_DIR}/${file_n} -O ${OUTPUT_DIR}/${file_n}.txt -H ${OUTPUT_DIR}/${file_n}.pdf
done

