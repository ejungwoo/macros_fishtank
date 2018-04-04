#!/usr/bin/env bash
#--- sbatch option ---#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --array=0-14

source /mnt/spirit/analysis/user/leej/SpiRITROOT.develop/build/config.sh
cd /mnt/spirit/analysis/user/leej/macros/reconstruction/

RUN=2900
NTOTAL=86916
NSPLIT=2000
GCData=/mnt/spirit/rawdata/misc/gainCalibration_groundPlane_120fC_117ns_20160509.root
GGData=/mnt/spirit/rawdata/misc/ggNoise/ggNoise_2895.root

SPLIT=$((3*SLURM_ARRAY_TASK_ID+0)); root run_reco_experiment_fieldMap.C\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\"$GCData\",\"$GGData\"\) -b -q -l > log/log_run$RUN\_$SPLIT.log 2>&1 &
SPLIT=$((3*SLURM_ARRAY_TASK_ID+1)); root run_reco_experiment_fieldMap.C\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\"$GCData\",\"$GGData\"\) -b -q -l > log/log_run$RUN\_$SPLIT.log 2>&1 &
SPLIT=$((3*SLURM_ARRAY_TASK_ID+2)); root run_reco_experiment_fieldMap.C\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\"$GCData\",\"$GGData\"\) -b -q -l > log/log_run$RUN\_$SPLIT.log 2>&1 &

wait
