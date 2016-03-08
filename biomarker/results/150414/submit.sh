#!/bin/bash

PREFIX="jjk_pipeline_test_Cancer"
VAR=$RANDOM
NCV=10

tDate=`date +%y%m%d`
tDir='/protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/results/'
tDir="$tDir$tDate/$PREFIX"
mkdir -p $tDir
mkdir $tDir/logs
mkdir $tDir/snapshots
mkdir $tDir/cv
mkdir $tDir/fs
mkdir $tDir/summary
mkdir $tDir/concat

cp /protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/results/150414/nonmodular_pipeline.R $tDir/snapshots/$PREFIX"_pipeline_snapshot.R"

for i in `seq 1 $NCV`
do
  ID=$PREFIX"_fold_"$i"_"$VAR
  LOG="$tDir/logs/$ID.qlog"
  qsub -N $ID -o $LOG /protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/results/150414/submit.qsub $i $PREFIX $tDir
done
