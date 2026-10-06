#!/usr/bin/env bash
#SBATCH --job-name=refbasedassembly                   # Job name
#SBATCH --output=refbasedassembly.%a.%A.out           # File to which stdout will be written
#SBATCH --error=refbasedassembly.%a.%A.err            # File to which stderr will be written
#SBATCH --partition=court                             # Partition
#SBATCH --cpus-per-task=16                            # Number of cpus ask per task
#SBATCH --time=02-00:00                               # Runtime in DD-HH:MM
#SBATCH --mem=200G                                    # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                               # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr           # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# author: Nicolas Nesi, Lokdeep Teekas
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 08/07/2026
# ---------------------------------

# ---------------------------------
# Submission:
# sbatch --array=1-x%4 P01_ReferenceBasedAssembly.sh -g Rferrumequinum
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# conda switch environments
conda_switch() {
    set +u
    conda deactivate 2>/dev/null || true
    conda activate "$1"
    set -u
}

# Positional parameters
GROUP="ViroCaen"
PROJECT="SugarBats"
QSCORE="15"
GENOME_REF=""

usage() {
  local exit_code="${1:-1}"
  {
    echo "Usage: sbatch --array=1-x $0 -g GENOME_REF" >&2
    echo "  -g|--genomeref Mandatory. Name of the genome of reference for aligment" >&2
    echo "  -h|--help  Show this help"
   } >&2 
    exit "$exit_code"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -g|--genomeref) GENOME_REF="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    *) 
         echo "Invalid option: $1" >&2
         usage 1
         ;;
    esac
done

if [[ -z "$GENOME_REF" ]]; then
    echo "Erreur: -g is mandatory" >&2
    exit 1
fi

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
BUSCO_DB="/dlocal/home/2019013/Databases/BUSCO/lineages/mammalia_odb12"
GENOME_DIR="/dlocal/home/2019013/Databases/bwamem_index/Mammals"
READS="$FOLDER/InputIlluminaReads"
OUTPUTREADS="$FOLDER/OutputCleanedReads"
OUTPUTASSEMBLY="$FOLDER/OutputAssembly"
OUTPUTBUSCO="$FOLDER/OutputBUSCO"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "$FOLDER/Scripts/List_samples.txt" | awk '{print $1}')"
SPNAME="$(grep -i "^$SAMPLE\b" "$FOLDER/Scripts/List_samples.txt" | awk '{print $2}')"

# Variables
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"

# Echo variables
echo "sample $SAMPLE"
echo "species name $SPNAME"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${READS}/${SPNAME}_${SAMPLE}"_R*.fastq.gz "${LOCAL_WORK_DIR}"

cd "${LOCAL_WORK_DIR}"

gunzip *.fastq.gz

echo "Analysis of ${SPNAME} genome: ${SAMPLE} run the ${RUNDATE} in Working directory ${PWD} job ID is ${SLURM_JOB_ID} and task of the array is ${SLURM_ARRAY_TASK_ID}"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
conda activate mapping_env
set -u
# include:
# adapterremoval=3.0.1
# bcftools=1.24
# bwamem2=2.3
# fastp=1.3.6
# samtools=1.24
# ---------------------------------

# ---------------------------------
# Raw reads cleaning
# ---------------------------------
echo "Starting raw reads cleaning"

adapterremoval3 \
--in-file1 "${SPNAME}_${SAMPLE}_R1.fastq" \
--in-file2 "${SPNAME}_${SAMPLE}_R2.fastq" \
--threads "${THREADS}" \
--out-json "${SPNAME}_${SAMPLE}"_Report.json \
--out-html "${SPNAME}_${SAMPLE}"_Report.html \
--out-file1 "${SPNAME}_${SAMPLE}_CleanedReads_R1.fastq.gz" \
--out-file2 "${SPNAME}_${SAMPLE}_CleanedReads_R2.fastq.gz"
# ---------------------------------

# ---------------------------------
# Reference based mapping
# ---------------------------------
echo "Starting mapping"

if [[ ! -f "${GENOME_DIR}/${GENOME_REF}.bwt.2bit.64" ]]; then
   bwa-mem2 index -p "${GENOME_DIR}/${GENOME_REF}" "${GENOME_DIR}/${GENOME_REF}"
   samtools faidx "${GENOME_DIR}/${GENOME_REF}"
fi

# Add a read group (helps downstream consistency; harmless if single sample)
RG="@RG\tID:${SPNAME}_${SAMPLE}\tSM:${SPNAME}_${SAMPLE}\tPL:ILLUMINA"

bwa-mem2 mem \
-t "${THREADS}" \
-R "${RG}" \
"${GENOME_DIR}/${GENOME_REF}" \
"${SPNAME}_${SAMPLE}_CleanedReads_R1.fastq.gz" "${SPNAME}_${SAMPLE}_CleanedReads_R2.fastq.gz" \
> "${SPNAME}_${SAMPLE}_${GENOME_REF}.sam"
# ---------------------------------

# ---------------------------------
# BAM processing
# ---------------------------------
echo "Starting BAM processing"

samtools view -@ "${THREADS}" -bS "${SPNAME}_${SAMPLE}_${GENOME_REF}.sam" > "${SPNAME}_${SAMPLE}_${GENOME_REF}.unsorted.bam"
rm -f "${SPNAME}_${SAMPLE}_${GENOME_REF}.sam"

samtools collate -@ "${THREADS}" -o "${SPNAME}_${SAMPLE}_${GENOME_REF}.collated.bam" "${SPNAME}_${SAMPLE}_${GENOME_REF}.unsorted.bam"
rm -f "${SPNAME}_${SAMPLE}_${GENOME_REF}.unsorted.bam"

samtools fixmate -@ "${THREADS}" -m "${SPNAME}_${SAMPLE}_${GENOME_REF}.collated.bam" "${SPNAME}_${SAMPLE}_${GENOME_REF}.fixmated.bam"
rm -f "${SPNAME}_${SAMPLE}_${GENOME_REF}.collated.bam"

samtools sort -@ "${THREADS}" -m 4G -o "${SPNAME}_${SAMPLE}_${GENOME_REF}.positionsort.bam" "${SPNAME}_${SAMPLE}_${GENOME_REF}.fixmated.bam"
rm -f "${SPNAME}_${SAMPLE}_${GENOME_REF}.fixmated.bam"

samtools markdup -@ "${THREADS}" -s "${SPNAME}_${SAMPLE}_${GENOME_REF}.positionsort.bam" "${SPNAME}_${SAMPLE}_${GENOME_REF}.mkdup.bam"
rm -f "${SPNAME}_${SAMPLE}_${GENOME_REF}.positionsort.bam"

samtools index -@ "${THREADS}" -c "${SPNAME}_${SAMPLE}_${GENOME_REF}.mkdup.bam"
# ---------------------------------

# ---------------------------------
# Variant calling
# ---------------------------------
# Strategy:
# - Call variants normally
# - Normalize
# - Instead of DROPPING lots of sites, SET low-confidence genotypes to missing (.)
#   This prevents tiny VCFs and avoids forcing reference where uncertain in consensus (since -M N uses N for missing).

bcftools mpileup \
--threads "${THREADS}" \
-Ou \
-f "${GENOME_DIR}/${GENOME_REF}" \
--min-MQ 20 --min-BQ 15 -d 10000 \
-a DP,AD,ADF,ADR \
"${SPNAME}_${SAMPLE}_${GENOME_REF}.mkdup.bam" \
| bcftools call -m -Ou --ploidy 2 \
| bcftools norm -Ou -f "${GENOME_DIR}/${GENOME_REF}" -m -any \
| bcftools +setGT -Ou -- \
-t q \
-n . \
-i 'QUAL<20 || FMT/DP<6' \
| bcftools view -Oz -o "${SPNAME}_${SAMPLE}.filtGT.vcf.gz"

bcftools index -f "${SPNAME}_${SAMPLE}.filtGT.vcf.gz"
# ---------------------------------

# ---------------------------------
# Consensus
# ---------------------------------
bcftools consensus \
-f "${GENOME_DIR}/${GENOME_REF}" \
-M N \
--iupac-codes \
"${SPNAME}_${SAMPLE}.filtGT.vcf.gz" \
> "${SPNAME}_${SAMPLE}_consensus.fa"

## Correcting weird symbols in the genome
awk '/^>/ { print; next } { gsub(/[^ACGTNacgtn]/,"N"); print $0 }' "${SPNAME}_${SAMPLE}_consensus.fa" > temp.fa
mv temp.fa "${SPNAME}_${SAMPLE}_consensus.fa"
# ---------------------------------

# ---------------------------------
# Re mapping to newly generated consensus
# ---------------------------------
GENOME_DIR="${PWD}"
NEWGENOME=$(ls | grep "${SPNAME}_${SAMPLE}_consensus.fa" | head -n 1)

# Add a read group (helps downstream consistency; harmless if single sample)
RG="@RG\tID:${SPNAME}_${SAMPLE}\tSM:${SPNAME}_${SAMPLE}\tPL:ILLUMINA"

bwa-mem2 index "${NEWGENOME}"

bwa-mem2 mem \
-t "${THREADS}" \
-R "$RG" \
"${GENOME_DIR}/${NEWGENOME}" \
"${SPNAME}_${SAMPLE}_CleanedReads_R1.fastq.gz" "${SPNAME}_${SAMPLE}_CleanedReads_R2.fastq.gz" \
> "${SPNAME}_${SAMPLE}_${NEWGENOME}.sam"
# ---------------------------------

# ---------------------------------
# BAM processing
# ---------------------------------
## Indexing for samtools
samtools faidx "${NEWGENOME}"
 
samtools view -@ "${THREADS}" -bS "${SPNAME}_${SAMPLE}_${NEWGENOME}.sam" > "${SPNAME}_${SAMPLE}_${NEWGENOME}.unsorted.bam"
rm -f "${SPNAME}_${SAMPLE}_${NEWGENOME}.sam"
 
samtools collate -@ "${THREADS}" -o "${SPNAME}_${SAMPLE}_${NEWGENOME}.collated.bam" "${SPNAME}_${SAMPLE}_${NEWGENOME}.unsorted.bam"
rm -f "${SPNAME}_${SAMPLE}_${NEWGENOME}.unsorted.bam"
 
samtools fixmate -@ "${THREADS}" -m "${SPNAME}_${SAMPLE}_${NEWGENOME}.collated.bam" "${SPNAME}_${SAMPLE}_${NEWGENOME}.fixmated.bam"
rm -f "${SPNAME}_${SAMPLE}_${NEWGENOME}.collated.bam"
 
samtools sort -@ "${THREADS}" -m 4G -o "${SPNAME}_${SAMPLE}_${NEWGENOME}.positionsort.bam" "${SPNAME}_${SAMPLE}_${NEWGENOME}.fixmated.bam"
rm -f "${SPNAME}_${SAMPLE}_${NEWGENOME}.fixmated.bam"
 
samtools markdup -@ "${THREADS}" -s "${SPNAME}_${SAMPLE}_${NEWGENOME}.positionsort.bam" "${SPNAME}_${SAMPLE}_${NEWGENOME}.mkdup.bam"
rm -f "${SPNAME}_${SAMPLE}_${NEWGENOME}.positionsort.bam"
 
samtools index -@ "${THREADS}" -c "${SPNAME}_${SAMPLE}_${NEWGENOME}.mkdup.bam"
 
## QC
samtools flagstat "${SPNAME}_${SAMPLE}_${NEWGENOME}.mkdup.bam" > "${SPNAME}_${SAMPLE}_${NEWGENOME}_mapping_stats.txt"
samtools stats "${SPNAME}_${SAMPLE}_${NEWGENOME}.mkdup.bam" > "${SPNAME}_${SAMPLE}_${NEWGENOME}.stats.txt"
# ---------------------------------

# ---------------------------------
# Variant calling
# ---------------------------------
# Strategy:
# - Call variants normally
# - Normalize
# - Instead of DROPPING lots of sites, SET low-confidence genotypes to missing (.)
#   This prevents tiny VCFs and avoids forcing reference where uncertain in consensus (since -M N uses N for missing).

bcftools mpileup \
--threads "${THREADS}" \
-Ou \
-f "${GENOME_DIR}/${NEWGENOME}" \
--min-MQ 20 --min-BQ 15 -d 10000 \
-a DP,AD,ADF,ADR \
"${SPNAME}_${SAMPLE}_${NEWGENOME}.mkdup.bam" \
| bcftools call -m -Ou --ploidy 2 \
| bcftools norm -Ou -f "${GENOME_DIR}/${NEWGENOME}" -m -any \
| bcftools +setGT -Ou -- \
-t q \
-n . \
-i 'QUAL<20 || FMT/DP<6' \
| bcftools view -Oz -o "${SPNAME}_${SAMPLE}_${NEWGENOME}.filtGT.vcf.gz"
 
bcftools index -f "${SPNAME}_${SAMPLE}_${NEWGENOME}.filtGT.vcf.gz"
# ---------------------------------
 
# ---------------------------------
# Consensus
# ---------------------------------
bcftools consensus \
-f "${GENOME_DIR}/${NEWGENOME}" \
-M N \
"${SPNAME}_${SAMPLE}_${NEWGENOME}.filtGT.vcf.gz" \
> "${SPNAME}_${SAMPLE}_${NEWGENOME}_remapped_consensus.fa"
 
## Correcting weird symbols in the genome
awk '/^>/ { print; next } { gsub(/[^ACGTNacgtn]/,"N"); print $0 }' "${SPNAME}_${SAMPLE}_${NEWGENOME}_remapped_consensus.fa" > temp.fa
mv temp.fa "${SPNAME}_${SAMPLE}_remapped_consensus.fa"
# ---------------------------------

# ---------------------------------
# BUSCO
# ---------------------------------
conda_switch busco_env
# include:
# busco=6.1.0
 
echo "Running BUSCO on ${SPNAME}_${SAMPLE}"

mkdir -p BUSCO_out/

busco \
--in "${SPNAME}_${SAMPLE}_remapped_consensus.fa" \
--out BUSCO_out/ \
--lineage_dataset "${BUSCO_DB}" \
--mode genome \
--cpu "${THREADS}" \
--force \
--offline \
--metaeuk \
--datasets_version odb12

busco --plot BUSCO_out/
# ---------------------------------

# ---------------------------------
# Move Outputs
# ---------------------------------
echo "Move cleaned reads"

mv "${LOCAL_WORK_DIR}"/"${SPNAME}_${SAMPLE}_CleanedReads_R1.fastq.gz" "${OUTPUTREADS}"

mv "${LOCAL_WORK_DIR}"/"${SPNAME}_${SAMPLE}_CleanedReads_R2.fastq.gz" "${OUTPUTREADS}"

mv "${SPNAME}_${SAMPLE}_Report.json" "${OUTPUTREADS}"

mv "${SPNAME}_${SAMPLE}_Report.html" "${OUTPUTREADS}"

echo "Move consensus assembly"

mv "${LOCAL_WORK_DIR}"/"${SPNAME}_${SAMPLE}_remapped_consensus.fa" "${OUTPUTASSEMBLY}/${SPNAME}_${SAMPLE}_final.fa"

pigz -p "${THREADS}" "${OUTPUTASSEMBLY}/${SPNAME}_${SAMPLE}_final.fa"

mv "${LOCAL_WORK_DIR}"/"${SPNAME}_${SAMPLE}_${NEWGENOME}.filtGT.vcf.gz" "${OUTPUTASSEMBLY}"

mv "${LOCAL_WORK_DIR}"/"${SPNAME}_${SAMPLE}_${NEWGENOME}_mapping_stats.txt" "${OUTPUTASSEMBLY}"

mv "${LOCAL_WORK_DIR}"/"${SPNAME}_${SAMPLE}_${NEWGENOME}.stats.txt" "${OUTPUTASSEMBLY}"

echo "Move BUSCO results"

mv "${LOCAL_WORK_DIR}/BUSCO_out/run_mammalia_odb12/short_summary.txt" "${OUTPUTBUSCO}/${SPNAME}_${SAMPLE}_short_summary.txt"

mv "${LOCAL_WORK_DIR}/BUSCO_out/run_mammalia_odb12/full_table.tsv" "${OUTPUTBUSCO}/${SPNAME}_${SAMPLE}_full_table.tsv"

mv "${LOCAL_WORK_DIR}/BUSCO_out/busco_figure.png" "${OUTPUTBUSCO}/${SPNAME}_${SAMPLE}_busco_figure.png"

cp -r "${LOCAL_WORK_DIR}/BUSCO_out/run_mammalia_odb12/busco_sequences/single_copy_busco_sequences/" "${OUTPUTBUSCO}/${SPNAME}_${SAMPLE}_single_copy_busco_sequences/"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/refbasedassembly."${SLURM_ARRAY_TASK_ID}".*.err "${FOLDER}/Scripts/Sbatch_log"
# ---------------------------------
