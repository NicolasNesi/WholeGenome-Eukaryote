#!/usr/bin/env bash
#SBATCH --job-name=LASTZ                              # Job name
#SBATCH --output=LASTZ.%a.%A.out                      # File to which stdout will be written
#SBATCH --error=LASTZ.%a.%A.err                       # File to which stderr will be written
#SBATCH --partition=long                              # Partition
#SBATCH --cpus-per-task=16                            # Number of cpus ask per task
#SBATCH --time=03-00:00                               # Runtime in DD-HH:MM
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
# sbatch --array=1-x P03_LASTZ.sh
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
LASTZ="/soft/2019013/Logiciels/make_lastz_chains"
REFGENOME="Homo_sapiens.GRCh38.filtered"
INPUTMASKING="$FOLDER/OutputMaskingGenome"
OUTPUTLASTZ="$FOLDER/OutputLastz"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "$FOLDER/Scripts/List_samples.txt" | awk '{print $1}')"
SPNAME="$(grep -i "^$SAMPLE\b" "$FOLDER/Scripts/List_samples.txt" | awk '{print $2}')"

# Variables
shortnameref="$(echo $REFGENOME|sed 's/.GRCh38.filtered//g')"
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
module load tools/nextflow/25.10.4
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${INPUTMASKING}/${SPNAME}_${shortnamesample}/${SPNAME}_${shortnamesample}.filtered.fa.masked.wm" "${LOCAL_WORK_DIR}"

cp "/dlocal/home/2019013/Databases/Ensembl_Genomes/Homo_sapiens.GRCh38.filtered.2bit" "${LOCAL_WORK_DIR}"

cp "/soft/2019013/Logiciels/make_lastz_chains/params.json" "${LOCAL_WORK_DIR}"

cd "${LOCAL_WORK_DIR}"

echo "Analysis of ${SPNAME} genome: sample ${SAMPLE} run the ${RUNDATE} in Working directory ${PWD} job ID is ${SLURM_JOB_ID} and task of the array is ${SLURM_ARRAY_TASK_ID}"
# ---------------------------------

# ---------------------------------
# Edit params.json
# ---------------------------------
sed \
  -e "s|\"reference_name\": null,|\"reference_name\": \"${shortnameref}\",|" \
  -e "s|\"query_name\": null,|\"query_name\": \"${SPNAME}_${shortnamesample}\",|" \
  -e "s|\"reference_genome\": null,|\"reference_genome\": \"${LOCAL_WORK_DIR}/${REFGENOME}.2bit\",|" \
  -e "s|\"query_genome\": null,|\"query_genome\": \"${LOCAL_WORK_DIR}/${SPNAME}_${shortnamesample}.filtered.fa.masked.wm\",|" \
  params.json > "${SPNAME}_params.json"
  # ---------------------------------

# ---------------------------------
# LASTZ
# ---------------------------------
echo "Check files in ${LOCAL_WORK_DIR}:"

ls -lh "${REFGENOME}.2bit" "${SPNAME}_${shortnamesample}.filtered.fa.masked.wm"

echo "Starting genome pariwise alignment with lastz"

mkdir -p results/

nextflow run /soft/2019013/Logiciels/make_lastz_chains/main.nf \
-params-file "${SPNAME}_params.json" \
-profile slurm,apptainer
# ---------------------------------

# ---------------------------------
# Move Outputs
# ---------------------------------
echo "Move results files"

# Create target directory
mkdir -p "${OUTPUTLASTZ}/${SPNAME}_${shortnamesample}/"

# Move and rename files with checks
echo "Moving *final.chain* files..."

ls "${LOCAL_WORK_DIR}"/results/07_final/*

file=$(ls "${LOCAL_WORK_DIR}"/results/07_final/*chain.gz 2>/dev/null)
if [ -f "$file" ]; then
   mv "$file" "${OUTPUTLASTZ}/${SPNAME}_${shortnamesample}/"
else 
   echo "No final.chain.gz generated !"
fi

echo "Moving 2bit files"

file="${LOCAL_WORK_DIR}/results/00_genome_prep/${SPNAME}_${shortnamesample}.2bit"
if [ -f "$file" ]; then
    mv "$file" "${INPUTMASKING}/${SPNAME}_${shortnamesample}/"
else
    echo "No 2bit file generated !"    
fi

echo "Moving reports"

if compgen -G "results/pipeline_info/"* > /dev/null; then
    mv results/pipeline_info/* "${OUTPUTLASTZ}/${SPNAME}_${shortnamesample}/"
else
    echo "Warning: report not found."
fi

echo "Moving nextflow.log and params.json"

mv .nextflow.log "${OUTPUTLASTZ}/${SPNAME}_${shortnamesample}/${SPNAME}_${shortnamesample}.nextflow.log"

cp "${SPNAME}_params.json" "${OUTPUTLASTZ}/${SPNAME}_${shortnamesample}"

echo "Finish!"
# ---------------------------------

# ---------------------------------
# Chaintools stats
# ---------------------------------
module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
conda activate chains_env
set -u
# include:
# chaintools=0.0.12

cd "${OUTPUTLASTZ}/${SPNAME}_${shortnamesample}/"

chaintools stats \
--threads "${THREADS}" \
--chain *.allfilled.chain.gz > "${SPNAME}_${shortnamesample}_chains_stats.tsv"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/LASTZ."${SLURM_ARRAY_TASK_ID}".*.err "${FOLDER}/Scripts/Sbatch_log"
# ---------------------------------
