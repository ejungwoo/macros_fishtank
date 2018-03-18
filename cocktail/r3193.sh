#!/usr/bin/env bash
#--- sbatch option ---#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-2

source /mnt/spirit/analysis/user/leej/SpiRITROOT.develop/build/config.sh
cd /mnt/spirit/analysis/user/leej/macros/cocktail/

RUN=3193
NTOTAL=78555
NSPLIT=10000
GCData=
GGData=

SPLIT=$((3*SLURM_ARRAY_TASK_ID+0)); root run_reco_experiment_fieldMap.C\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\"$GCData\",\"$GGData\"\) -b -q -l > data/log_run$RUN\_$SPLIT.log 2>&1 &
SPLIT=$((3*SLURM_ARRAY_TASK_ID+1)); root run_reco_experiment_fieldMap.C\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\"$GCData\",\"$GGData\"\) -b -q -l > data/log_run$RUN\_$SPLIT.log 2>&1 &
SPLIT=$((3*SLURM_ARRAY_TASK_ID+2)); root run_reco_experiment_fieldMap.C\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\"$GCData\",\"$GGData\"\) -b -q -l > data/log_run$RUN\_$SPLIT.log 2>&1 &

wait
