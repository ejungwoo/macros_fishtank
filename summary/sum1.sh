#!/usr/bin/env bash
#--- sbatch option ---#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000

source /mnt/spirit/analysis/user/leej/SpiRITROOT.develop/build/config.sh
cd /mnt/spirit/analysis/user/leej/macros/summary/

root -b -q make_summary.C\(1\) > data/batch1.log
