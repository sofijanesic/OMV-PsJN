#!/bin/sh
#
# Run ....
#
#SBATCH -J umiT-fastp                     # A single job name for the array
#SBATCH --cpus-per-task 8           # 8 CPUs (threads)
#SBATCH --nodes 1                    # one node
#SBATCH --time 2-00                  # Running time of 2 days
#SBATCH --mem 32G                    # Memory request
#SBATCH --output /lustre/imgge/grupa08/OMV/logs/%x_%A_%a.out       # Standard output
#SBATCH --error  /lustre/imgge/grupa08/OMV/logs/%x_%A_%a.err       # Standard error
#SBATCH --array=0-7%4                  # radi 4 odjednom

module purge
module load bioinfo
module load fastp
module load UMI-tools

SAMPLES="CELLS1_SE50 CELLS3_SE50 OMV1_SE50 OMV3_SE50 CELLS1_PE100 CELLS3_PE100 OMV1_PE100 OMV3_PE100"

set -ex

# split samples into array
IFS=' ' read -r -a samples_arr <<< "$SAMPLES"

SAMPLE="${samples_arr[$SLURM_ARRAY_TASK_ID]}"


BASE="/lustre/imgge/grupa08/OMV"
UMI_EXTRACTION=false
DATA="${BASE}/fastq"
EXTRACTED="${BASE}/extracted"
TRIMMED="${BASE}/trimmed"


if UMI_EXTRACTION; then
# UMIs EXTRACTION
umi_tools extract \
	--ignore-read-pair-suffixes \
	--extract-method=regex \
	--bc-pattern2='.*(?P<umi_1>.{10})' \
	--read2-only \
	--stdin="${DATA}/${SAMPLE}_1.fq.gz" \
	--read2-in="${DATA}/${SAMPLE}_2.fq.gz" \
	--stdout="${EXTRACTED}/${SAMPLE}_extracted_1.fq.gz" \
	--read2-out="${EXTRACTED}/${SAMPLE}_extracted_2.fq.gz" \
	--log2stderr
fi

# TRIMMING
quality_opts="--overrepresentation_analysis --qualified_quality_phred 20 --unqualified_percent_limit 30 --average_qual 25"
length_opts="--length_required 15"
adapters_opt="--adapter_sequence=AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
cut_right_opts="--cut_right_window_size 4 --cut_right_mean_quality 15 --cut_right"

r2_files=""
if [[ "$SAMPLE" == *"PE100"* ]]; then
	r2_files="-I ${EXTRACTED}/${SAMPLE}_extracted_2.fq.gz -O ${TRIMMED}/${SAMPLE}_trimmed_2.fq.gz"
	length_opts="--length_required 15"
fi

fastp -w $SLURM_CPUS_PER_TASK $quality_opts $length_opts $cut_right_opts \
	-i "${EXTRACTED}/${SAMPLE}_extracted_1.fq.gz" -o "${TRIMMED}/${SAMPLE}_trimmed_1.fq.gz" \
	$r2_files \
	-j "${TRIMMED}/${SAMPLE}_fastp.json" -h "${TRIMMED}/$reports_dir/${SAMPLE}_fastp.html" \
	$adapters_opt
