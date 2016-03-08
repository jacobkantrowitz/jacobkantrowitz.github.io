#!/bin/bash

PREFIX="jjk_pipeline_test_Cancer"
VAR=$RANDOM
NCV=10

tDate=`date +%y%m%d`
tDir='/protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/results/'
tDir="$tDir$tDate/$PREFIX"


PPLOG="$tDir/logs/"$PREFIX"_post_process.log"

qsub -N "postProcess" -o $PPLOG -hold_jid $PREFIX"_fold_*" /protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/results/150414/postProcess.qsub $NCV $PREFIX $tDir





