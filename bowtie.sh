THREADS=4   # use half of 8 cores
REF=GCA_000020125.1_ASM2012v1_genomic.fna

# Step 1: Build Bowtie index (only once per reference)
bowtie-build $REF reference_index

# Step 2: Align reads with Bowtie
bowtie -p $THREADS \
  reference_index \
  OMV1_PE100_1.fq.gz \
  OMV1_PE100_1_bowtie.sam

# Step 3: Convert SAM to BAM, sort, index, and run flagstat
samtools view -bS OMV1_PE100_1_bowtie.sam | \
samtools sort -o OMV1_PE100_1_bowtie_sorted.bam

samtools index OMV1_PE100_1_bowtie_sorted.bam
samtools flagstat OMV1_PE100_1_bowtie_sorted.bam
