#!/usr/bin/env bash
#SBATCH --job-name=RepeatMasking                       # Job name
#SBATCH --output=RepeatMasking.%a.%A.out               # File to which stdout will be written
#SBATCH --error=RepeatMasking.%a.%A.err                # File to which stderr will be written
#SBATCH --partition=long                               # Partition
#SBATCH --cpus-per-task=24                             # Number of cpus ask per task
#SBATCH --time=04-04:00                                # Runtime in DD-HH:MM
#SBATCH --mem=200G                                     # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                                # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr            # Email address

# ---------------------------------
# Version v1.1
# ---------------------------------
# author: Nicolas Nesi, Lokdeep Teekas
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 12/08/2026
# ---------------------------------

# ---------------------------------
# Submission:
# sbatch --array=1-x RepeatMasking.sh -f $FORMAT -t $TOP
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# Positional parameters
GROUP="ViroCaen"
PROJECT="SugarBats"
FORMAT=""
TOP=""

usage() {
    echo "Usage: $0 -f FORMAT -t TOP" >&2
	echo "   -f, --format  Mandatory. Either NCBI or bwamem" >&2
	echo "   -t, --top     Optional. Number of top (longest) scaffolds to keep." >&2
	echo "                 Omit this option entirely to keep all scaffold (do NOT pass -t with no value)." >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -f|--format) FORMAT="$2"; shift 2 ;;
	-t|--top) TOP="$2"; shift 2 ;;
    *) 
         echo "Invalid option: $1" >&2
         usage
         ;;
    esac
done

if [[ -z "$FORMAT" ]]; then
    echo "Erreur: -f is mandatory; either NCBI or bwamem" >&2
    exit 1
fi

# reject anything that isn't one of the two accepted values
if [[ "$FORMAT" != "NCBI" && "$FORMAT" != "bwamem" ]]; then
    echo "Erreur: -f must be NCBI or bwamem" >&2
    exit 1
fi

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATHSIF="/soft/2019013/Logiciels/TETools"
INPUTGENOME="$FOLDER/OutputAssembly"
OUTPUTMASKING="$FOLDER/OutputMaskingGenome"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "$FOLDER/Scripts/List_samples.txt" | awk '{print $1}')"
SPNAME="$(grep -i "^$SAMPLE\b" "$FOLDER/Scripts/List_samples.txt" | awk '{print $2}')"

# Variables
if [[ "$FORMAT" == "NCBI" ]]; then
   shortnamesample="$(echo "$SAMPLE" | sed 's/_genomic.fna//')"
elif [[ "$FORMAT" == "bwamem" ]]; then
   shortnamesample="$(echo "$SAMPLE")"
fi

RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"

# Echo variables
echo "sample $SAMPLE"
echo "species name $SPNAME"
echo "name sample $shortnamesample"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
conda activate NN_env
set -u
# include:
# python=3.13.1
# seqkit=2.13.0
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
if [[ "$FORMAT" == "NCBI" ]]; then
  cp "${INPUTGENOME}"/NCBI/"${SAMPLE}".gz "${LOCAL_WORK_DIR}"
elif [[ "$FORMAT" == "bwamem" ]]; then
  cp "${INPUTGENOME}"/"${SPNAME}_${SAMPLE}_final.fa.gz" "${LOCAL_WORK_DIR}"
fi

cd "${LOCAL_WORK_DIR}"

if [[ "$FORMAT" == "NCBI" ]]; then
 gunzip "${SAMPLE}".gz
elif [[ "$FORMAT" == "bwamem" ]]; then
 gunzip "${SPNAME}_${SAMPLE}_final.fa.gz" 
fi

echo "Analysis of ${SPNAME} genome : ${SAMPLE} run the ${RUNDATE} in Working directory ${PWD} job ID is ${SLURM_JOB_ID} and task of the array is ${SLURM_ARRAY_TASK_ID}"
# ---------------------------------

# ---------------------------------
# Cleaning and Masking of query genome
# ---------------------------------
echo "Filter scaffolds of $SAMPLE"

set +o pipefail

if [[ "$FORMAT" == "NCBI" && -n "$TOP" ]]; then
 awk '{print $1}' "${SAMPLE}" | seqkit sort -l -r | seqkit head -n "${TOP}" | sed -e 's/\./_/g' -e 's/|/_/g' -e 's/:/_/g' -e 's/-/_/g' > "${shortnamesample}.filtered.fa"
elif [[ "$FORMAT" == "NCBI" && -z "$TOP" ]]; then
 awk '{print $1}' "${SAMPLE}" | seqkit sort -l -r | sed -e 's/\./_/g' -e 's/|/_/g' -e 's/:/_/g' -e 's/-/_/g' > "${shortnamesample}.filtered.fa"
elif [[ "$FORMAT" == "bwamem" && -n "$TOP" ]]; then
 awk '{print $1}' "${SPNAME}_${SAMPLE}_final.fa" | seqkit sort -l -r | seqkit head -n "${TOP}" | sed -e 's/\./_/g' -e 's/|/_/g' -e 's/:/_/g' -e 's/-/_/g' > "${shortnamesample}.filtered.fa"
elif [[ "$FORMAT" == "bwamem" && -z "$TOP" ]]; then
 awk '{print $1}' "${SPNAME}_${SAMPLE}_final.fa" | seqkit sort -l -r | sed -e 's/\./_/g' -e 's/|/_/g' -e 's/:/_/g' -e 's/-/_/g' > "${shortnamesample}.filtered.fa"
fi

set -o pipefail
 
if [ ! -f "${shortnamesample}.filtered.fa" ]; then
    echo "Error: ${shortnamesample}.filtered.fa failed!"
    exit 1
fi

echo "Starting masking of $SAMPLE"

singularity run "${PATHSIF}"/dfam-tetools-latest.sif \
BuildDatabase \
-name "${shortnamesample}.filtered" \
"${shortnamesample}.filtered.fa"

singularity run "${PATHSIF}"/dfam-tetools-latest.sif \
RepeatModeler \
-database "${shortnamesample}.filtered" \
--threads "${THREADS}" \
-LTRStruct

singularity run "${PATHSIF}"/dfam-tetools-latest.sif \
RepeatMasker \
-pa "${THREADS}" \
-lib "${shortnamesample}.filtered-families.fa" \
-xsmall \
-q \
-gff \
"${shortnamesample}.filtered.fa"

singularity run "${PATHSIF}"/dfam-tetools-latest.sif \
windowmasker \
-mk_counts \
-in "${shortnamesample}.filtered.fa" \
-out "${shortnamesample}.filtered.wm.counts"

singularity run "${PATHSIF}"/dfam-tetools-latest.sif \
windowmasker \
-ustat "${shortnamesample}.filtered.wm.counts" \
-in "${shortnamesample}.filtered.fa.masked" \
-out "${shortnamesample}.filtered.fa.masked.wm" \
-outfmt fasta
# ---------------------------------

# ---------------------------------
# Move Outputs
# ---------------------------------
echo "Move outputs files"

# Create target directory
if [[ ! -d "${OUTPUTMASKING}/${SPNAME}_${shortnamesample}/" ]]; then
 mkdir -p "${OUTPUTMASKING}/${SPNAME}_${shortnamesample}/"
fi

# Move and rename files with checks
echo "Moving ${shortnamesample}.filtered.fa.masked.wm"

if [ -f "${shortnamesample}.filtered.fa.masked.wm" ]; then
    mv "${shortnamesample}.filtered.fa.masked.wm" "${OUTPUTMASKING}/${SPNAME}_${shortnamesample}/${SPNAME}_${shortnamesample}.filtered.fa.masked.wm"
else
    echo "Warning: ${shortnamesample}.filtered.fa.masked.wm"
fi

echo "Done."
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
#rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/RepeatMasking."${SLURM_ARRAY_TASK_ID}".*.err "${FOLDER}/Scripts/Sbatch_log"
# ---------------------------------

