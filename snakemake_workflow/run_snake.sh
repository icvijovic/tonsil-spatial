#!/bin/bash

#name shell variables for calling in snakemake and sbatch
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
RESTART=0

SNAKEFILE=Snakefile.smk

#Snakemake config
NJOBS=200
WAIT=120

if [ $USER = "mswift2" ]
    then
    source /home/groups/quake/mswift/mambaforge/etc/profile.d/conda.sh
    conda activate scanpy_latest
    mkdir -p snakemake_logs/slurm_logs/
elif [ $USER = "cvijovic" ]
    then
    source /home/groups/quake/cvijovic/miniconda3/etc/profile.d/conda.sh
    conda activate snakemake
    mkdir -p snakemake_logs/
else
    mkdir -p snakemake_logs/
    echo "unknown user; will not explicitly activate conda"
fi

#log file for process that calls snakemake
SBATCH_LOGFILE=snakemake_logs/cluster.$DATETIME.log
SBATCH_LOGFILE_ERR=snakemake_logs/cluster.$DATETIME.log.err

# for a specific cell / output
TARGET=''
if [ $# -eq 0 ]
  then
    # Dry run snakemake
    snakemake -s $SNAKEFILE $TARGET --use-conda --keep-target-files --rerun-incomplete -n -r --quiet --keep-going --rerun-triggers mtime

elif [ $1 = "unlock" ]
    then
        snakemake -s $SNAKEFILE $TARGET -F --rerun-incomplete --unlock --cores 1

elif [ $1 = "touch" ]
    then
        snakemake -s $SNAKEFILE $TARGET -F --rerun-incomplete --unlock --touch --cores 1
    
elif [ $1 = "dry" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --profile profile/ -n

elif [ $1 = "profile" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --profile profile/

elif [ $1 = "sbatch" ]
    # Run snakemake as an SBATCH job, good for long workflows
    then
    sbatch \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=8000 \
        --mail-user=$EMAIL \
        --time 4-0 \
        -p quake \
        -o $SBATCH_LOGFILE \
        -e $SBATCH_LOGFILE_ERR \
        run_snake.sh profile

else
    echo "wrong option"
fi
