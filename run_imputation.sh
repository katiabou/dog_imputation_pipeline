#!/usr/bin/env bash

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# get the command line arguments
args=("$@")

datetime=$(date +%Y-%m-%d-%H%M)

# make a unique log file
logfile="dog_imputation-${datetime}.log"

# print the server name and start time to the log file
echo "SERVER: $HOSTNAME" >>${logfile}
echo "DATE: ${datetime}" >>${logfile}
#echo "TARGETS: " "${args[@]}" >>$logfile

# load the conda environment
eval "$(conda shell.bash hook)"
conda activate dogs

# maximum number of concurrent FTP requests (prevents overloading the FTP server)
#MAX_FTP=10

if ! command -v free &>/dev/null; then
   #MacOS does not have the free command
   MAX_MEM=$(sysctl -a | awk '/^hw.memsize:/{print $2/(1024)^2}')
else
   #but linux does
   MAX_MEM=$(free -m | awk '/^Mem:/{print $2}')
fi

# hide 1 GB of RAM from the snakemake scheduler, to avoid exhausting total system RAM
MAX_MEM=$((MAX_MEM - 1024))

flags="--cores all "
flags+="--nolock "
flags+="--keep-going "
flags+="--printshellcmds "
flags+="--show-failed-logs "
flags+="--keep-incomplete "
flags+="--rerun-incomplete "
flags+="--reason "
flags+="--resources mem_mb=${MAX_MEM} "

snakemake ${flags} -- "${args[@]}" &>>${logfile}
#snakemake -np &>>${logfile}

echo "DONE!" >>${logfile}

