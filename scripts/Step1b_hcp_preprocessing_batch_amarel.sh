#!/bin/bash

# This script runs a bash function on a large number of subjects (executing the HCP Preprocessing Pipeliens) using the supercomputer queuing system

#modified by R Mill for CPRO2_learning preproc
#adapted for Amarel
#**NOTE: make sure you scp the data out of scratch/*output/ once you get the email saying the job is finished!

#if re-running aborted subs from other servers make sure to scp:
#1) the raw input to /home/${usr_ID}/projects/CPRO2_learning/data/rawdata/mriqc_BIDS
#2) the output from the pipelines to /scratch/${usr_ID}/CPRO2_learning_${subjNum}/output/${subjNum}/.


##*** Modify these variables:
#RM NOTE - only need this var for setting up paths; not necessary for running on lab partition
usr_ID="f_keanebp"
#Where step1/1a/1b and opts.shlib should be
scriptDir="/projects/${usr_ID}/NeuralMechanisms/docs/scripts/" # updated by bpk, 8/3/18
#path to preproc script (Step1)
#preproc_script=${scriptDir}/Step1_hcp_preprocessing_msmall_amarel.sh
preproc_script=${scriptDir}/Step1_hcp_preprocessing_v2.sh
#where batch scripts for each subject are written
subjBatchScriptDir="${scriptDir}/subjbatch/"
if [ ! -e $subjBatchScriptDir ]; then mkdir -p $subjBatchScriptDir; fi
jobNamePrefix="H"

listOfSubjects="sub-C05"
# Completed subjects:
#

##Full list of subject numbers (keep full list of subject numbers commented-out for future reference):

##Make and execute a batch script for each subject
for subjNum in $listOfSubjects
do

 	cd ${subjBatchScriptDir}

	batchFilename=${subjNum}_hcppreprocBatch.sh

	#RM - last two parms send email when script ends, important so that one can scp promptly
	#modified from nm3: time (3-day limit, set to one day currently), partition (think only main is accessible for free?), cpus-per-task (was 20 on nm3;seems like most main nodes have max 28 cpus as shown by scontrol show node command )
	echo "#!/bin/bash" > $batchFilename
	echo "#SBATCH --time=02-00:00:00" >> $batchFilename
	echo "#SBATCH --nodes=1" >> $batchFilename
	echo "#SBATCH --ntasks=1" >> $batchFilename
	echo "#SBATCH --partition=p_keanebp" >> $batchFilename #
	echo "#SBATCH --job-name=${jobNamePrefix}${subjNum}" >> $batchFilename
	echo "#SBATCH --output=slurm.${jobNamePrefix}${subjNum}.out" >> $batchFilename
	echo "#SBATCH --error=slurm.${jobNamePrefix}${subjNum}.err" >> $batchFilename
	echo "#SBATCH --cpus-per-task=4" >> $batchFilename
	echo "#SBATCH --mail-type=END" >> $batchFilename
	echo "#SBATCH --mail-user=brian.keane@rutgers.edu" >> $batchFilename

	echo "#Run the MATLAB command" >> $batchFilename
	echo "cd $scriptDir" >> $batchFilename


	#command to execute HCP preproc script - all modules!
	#echo "${preproc_script} --server='amarel' --preFS='true' --FS='true' --postFS='true' --fmriVol='true' --fmriSurf='true' --restFix='true' --msmAll='true' --dedriftResample='true' --subj='${subjNum}'" >> $batchFilename
	echo "${preproc_script} --server='amarel' --preFS='false' --FS='false' --postFS='false' --fmriVol='true' --fmriSurf='true' --restFix='false' --msmAll='false' --dedriftResample='true' --subj='${subjNum}'" >> $batchFilename
	#echo "${preproc_script} --server='amarel' --fmriVol_test='true' --fmriVol_test_run='5' --fmriSurf='true' --restFix='true' --msmAll='true' --dedriftResample='true' --subj='${subjNum}'" >> $batchFilename
	#echo "${preproc_script} --server='amarel' --fmriVol='true' --fmriSurf='true' --restFix='true' --msmAll='true' --dedriftResample='true' --subj='${subjNum}'" >> $batchFilename

	#Submit the job
	sbatch $batchFilename

done
