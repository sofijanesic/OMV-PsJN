#!/bin/sh
#
# Run ....
#
#SBATCH -J bwa-dedup                     # A single job name for the array
#SBATCH --cpus-per-task 16           # 8 CPUs (threads)
#SBATCH --nodes 1                    # one node
#SBATCH --time 2-00                  # Running time of 2 days
#SBATCH --mem 64G                    # Memory request
#SBATCH --output /lustre/imgge/grupa08/OMV/logs/%x_%A_%a.out       # Standard output
#SBATCH --error  /lustre/imgge/grupa08/OMV/logs/%x_%A_%a.err       # Standard error
#SBATCH --array=0-7%4                  # radi 4 odjednom

module purge
module load bioinfo
module load BWA
module load SAMtools
module load UMI-tools

SAMPLES="CELLS1_SE50 CELLS3_SE50 OMV1_SE50 OMV3_SE50 CELLS1_PE100 CELLS3_PE100 OMV1_PE100 OMV3_PE100"

set -ex

# split samples into array
IFS=' ' read -r -a samples_arr <<< "$SAMPLES"

SAMPLE="${samples_arr[$SLURM_ARRAY_TASK_ID]}"


BASE='/lustre/imgge/grupa08/OMV'

data_dir="${BASE}"
TRIMMED="${BASE}/trimmed"
REF_FA="${BASE}/ref/genome.fa"
BAMS="${BASE}/bams"
DEDUP="${BASE}/dedup"



if [[ "$SAMPLE" == *"SE50"* ]]; then
# BWA ALN (SE)

bwa aln \
  -B 5 -O 4 -E 1 -l 12 \
  -t "$SLURM_CPUS_PER_TASK" \
  "$REF_FA" "${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" > "${BAMS}/${SAMPLE}.sai"

bwa samse \
  "$REF_FA" \
  "${BAMS}/${SAMPLE}.sai" \
  "${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" | \
samtools view -bS - | \
samtools sort -o "${BAMS}/${SAMPLE}_bwa_sorted.bam"
rm "${BAMS}/${SAMPLE}.sai"

elif [[ "$SAMPLE" == *"PE100"* ]]; then
# BWA ALN (PE)

bwa aln -B 5 -O 4 -E 1 -l 12 -t "$SLURM_CPUS_PER_TASK" \
  "$REF_FA" "${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" > "${BAMS}/${SAMPLE}_R1.sai"

bwa aln -B 5 -O 4 -E 1 -l 12 -t "$SLURM_CPUS_PER_TASK" \
  "$REF_FA" "${TRIMMED}/${SAMPLE}_trimmed_2.fq.gz" > "${BAMS}/${SAMPLE}_R2.sai"

bwa sampe \
  "$REF_FA" \
  "${BAMS}/${SAMPLE}_R1.sai" \
  "${BAMS}/${SAMPLE}_R2.sai" \
  "${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" \
  "${TRIMMED}/${SAMPLE}_trimmed_2.fq.gz" | \
samtools view -bS - | \
samtools sort -o "${BAMS}/${SAMPLE}_bwa_sorted.bam"
rm "${BAMS}/${SAMPLE}_R1.sai"
rm "${BAMS}/${SAMPLE}_R2.sai"

fi

samtools index "${BAMS}/${SAMPLE}_bwa_sorted.bam"


exit
# UMI DEDUP

umi_tools dedup \
--stdin "${BAMS}/${SAMPLE}_bwa_sorted.bam" \
--stdout "${DEDUP}/${SAMPLE}_bwa_dedup.bam" \
--method unique

samtools index "${DEDUP}/${SAMPLE}_bwa_dedup.bam"
