#!/usr/bin/env bash
#SBATCH --job-name=TOGA2                              # Job name
#SBATCH --output=TOGA2.%a.%A.out                      # File to which stdout will be written
#SBATCH --error=TOGA2.%a.%A.err                       # File to which stderr will be written
#SBATCH --partition=tcourt                            # Partition
#SBATCH --cpus-per-task=24                            # Number of cpus ask per task
#SBATCH --time=01-00:00                               # Runtime in DD-HH:MM
#SBATCH --mem=150G                                    # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                               # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr           # Email address

# ---------------------------------
# Version v1.4
# ---------------------------------
# author: Nicolas Nesi, Lokdeep Teekas
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 11/08/2026
# ---------------------------------

# ---------------------------------
# Submission:
# sbatch --array=1-x P04_TOGA2.sh
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# Positional parameters
GROUP="ViroCaen"
PROJECT="SugarBats"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
TOGA2="/soft/2019013/Logiciels/TOGA2_v209"
INPUTMASKING="$FOLDER/OutputMaskingGenome"
INPUTLASTZ="$FOLDER/OutputLastz"
OUTPUTTOGA2="$FOLDER/OutputTOGA2"
REFGENOME="Homo_sapiens.GRCh38.filtered"
nf_config_path="$TOGA2/nextflow_configs"
query_file_path="$FOLDER/OutputLastz"
ref_seq_path="/dlocal/home/2019013/Databases/Ensembl_Genomes"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "$FOLDER/Scripts/List_samples.txt" | awk '{print $1}')"
SPNAME="$(grep -i "^$SAMPLE\b" "$FOLDER/Scripts/List_samples.txt" | awk '{print $2}')"

# Variables
shortnameref="$(echo $REFGENOME|sed 's/.GRCh38.filtered.2bit//g')"
shortnamesample="$(echo $SAMPLE|sed 's/_genomic.fna//g')"
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"

# Echo variables
echo "sample $SAMPLE"
echo "species name $SPNAME"
echo "name ref $shortnameref"
echo "name sample $shortnamesample"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0 tools/nextflow/25.10.4
set +u
eval "$(conda shell.bash hook)"
conda activate toga2_env
set -u
# include:
# lastz=1.04.52
# python=3.11.14
# toga2=2.0.9f
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${INPUTMASKING}/${SPNAME}_${shortnamesample}/${SPNAME}_${shortnamesample}.2bit" "${LOCAL_WORK_DIR}"

cp "${INPUTLASTZ}/${SPNAME}_${shortnamesample}"/*.allfilled.chain.gz "${LOCAL_WORK_DIR}"

cp "/dlocal/home/2019013/Databases/Ensembl_Genomes/${REFGENOME}.2bit" "${LOCAL_WORK_DIR}"

cd "${LOCAL_WORK_DIR}"

query_chain="$(ls *.allfilled.chain.gz)"

echo "chain: $query_chain"

echo "Analysis of ${SPNAME} genome: ${SAMPLE} run the ${RUNDATE} in Working directory ${PWD} job ID is ${SLURM_JOB_ID} and task of the array is ${SLURM_ARRAY_TASK_ID}"
# ---------------------------------

# ---------------------------------
# TOGA2
# ---------------------------------
echo "Check files in ${LOCAL_WORK_DIR}:"

ls -lh "${REFGENOME}.2bit" "${SPNAME}_${shortnamesample}.2bit"

echo "Starting TOGA2"

mkdir -p output/

export TMPDIR="/soft/2019013/Logiciels/TOGA2_v209/tmp_toga2/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

python3 "${TOGA2}"/toga2.py run \
--ref_2bit "${REFGENOME}.2bit" \
--query_2bit "${SPNAME}_${shortnamesample}.2bit" \
--chain_file "${query_chain}" \
--ref_annotation "${ref_seq_path}/Homo_sapiens.GRCh38.bed" \
--nextflow_config_dir "${nf_config_path}" \
--project_name "${SPNAME}" \
--max_parallel_time 1 \
--no_isoform_file \
--no_u12_file \
--no_spliceai \
--ignore_crashed_parallel_batches \
--parallel_strategy slurm \
--output output/

grep -A11 "^Orthology resolution:" "${LOCAL_WORK_DIR}/output/summary.txt"
# ---------------------------------

# ---------------------------------
# Move Outputs
# ---------------------------------
echo "Move outputs files"

# Create target directory
mkdir -p "${OUTPUTTOGA2}/${SPNAME}_${shortnamesample}/"

cp -R "${LOCAL_WORK_DIR}"/output/* "${OUTPUTTOGA2}/${SPNAME}_${shortnamesample}/"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
rm -rf "$TMPDIR"
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/TOGA2."${SLURM_ARRAY_TASK_ID}".*.err "${FOLDER}/Scripts/Sbatch_log"
# ---------------------------------