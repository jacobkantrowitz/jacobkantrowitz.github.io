#!/bin/bash

# define the results directory, where items will be moved
results="/protected/projects/pulmarray/Allegro/COPD_Cancer/results/"

# define the experiments directory, where items will be identified and moved from
experiments="/protected/projects/pulmarray/Allegro/COPD_Cancer/experiments"

for file in `find $experiments -maxdepth 2 -not -name "*.R*" -and -not -name "*.R" -and -not -name "*.Rhistory" -and -not -name "*.qsub" -and -not -name "*.sh" -type f`
do
	# find last modified date of the file
	lastMod=`date -r $file +%Y-%m-%d`
	expmtDir=`dirname $file`
	expmtName=`basename $expmtDir`
	expmtResults=$results$expmtName/$lastMod
	#echo $expmtResults
	mkdir -p -v $expmtResults
	# create results directory with the given date if it doesn't already exist
	#dateRes=$results$lastMod
	#mkdir -p $dateRes
	# get experiment directory name
	#makeDir=`dirname $file | rev | cut -d'/' -f1 | rev`
	# within the specified date results directory, create a directory for the experiment
	#mkdir -p $dateRes/$makeDir
	# move all the files to the date/experiment directory
	mv $file $expmtResults

done




# define the name of the date-specific results directory to create if it doesn't exist
#dateRes=$results`date +%Y-%m-%d`
#mkdir -p $dateRes

#for file in `find $PWD -daystart -ctime 0 ! -name "*.Rmd" -type f`
#do
#	makeDir=`dirname $file | rev | cut -d'/' -f1 | rev`
#	echo "TEST"
#	mkdir -p $dateRes/$makeDir
#	mv $file $dateRes/$makeDir

#done

# do date -r $file +%Y-%m-%d

# for each file in experiments that is txt, pdf, etc

# create results directory with the given date if it doesn't already exist
# get experiment directory name
# within the specified date results directory, create a directory for the experiment
# move all the files to the date/experiment directory