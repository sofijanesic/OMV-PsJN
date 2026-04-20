#!/bin/sh
#
# Run ....
#
#SBATCH -J bbmap-dedup               # A single job name for the array
#SBATCH --cpus-per-task 16           # 8 CPUs (threads)
#SBATCH --nodes 1                    # one node
#SBATCH --time 2-00                  # Running time of 2 days
#SBATCH --mem 64G                    # Memory request
#SBATCH --output /lustre/imgge/grupa08/OMV/logs/%x_%A_%a.out       # Standard output
#SBATCH --error  /lustre/imgge/grupa08/OMV/logs/%x_%A_%a.err       # Standard error
#SBATCH --array=0-7%4                  # radi 4 odjednom

module purge
module load bioinfo
module load bbmap
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
# 5. BBMAP (SE50)

bbmap.sh \
  in="${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" \
  out="${BAMS}/${SAMPLE}_bb.sam" \
  ref="$REF_FA" \
  k=12 \
  minid=0.90 \
  maxindel=10 \
  threads=$SLURM_CPUS_PER_TASK

elif [[ "$SAMPLE" == *"PE100"* ]]; then
# === 5. BBMap mapping ===

bbmap.sh \
  in1="${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" \
  in2="${TRIMMED}/${SAMPLE}_trimmed_2.fq.gz" \
  out="${BAMS}/${SAMPLE}_bb.sam" \
  ref="$REF_FA" \
  k=12 \
  minid=0.90 \
  maxindel=10 \
  threads=$SLURM_CPUS_PER_TASK

fi


samtools view -bS "${BAMS}/${SAMPLE}_bb.sam" | \
samtools sort -o "${BAMS}/${SAMPLE}_bb_sorted.bam"
samtools index "${BAMS}/${SAMPLE}_bb_sorted.bam"

rm "${BAMS}/${SAMPLE}_bb.sam"

# UMI DEDUP

umi_tools dedup \
--stdin "${BAMS}/${SAMPLE}_bb_sorted.bam" \
--stdout "${DEDUP}/${SAMPLE}_bb_dedup.bam" \
--method unique

samtools index "${DEDUP}/${SAMPLE}_bb_dedup.bam"


