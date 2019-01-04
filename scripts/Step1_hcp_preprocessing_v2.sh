#!/bin/bash


#modified by Ravi Mill for CPRO2_learning
#assumes data are in BIDS format
# 03/16/18

# Modified again by Brian Keane to preprocess data from "Neural Mechanisms" study.
# 07-08/2018
# Changed anatomical scans to reconstruction 1 (rec-1)
# Changed field maps (SEPhaseNeg, SEPhasePos) to run-1 for Viscomp.

#**Note this version can only be used on Amarel (the directory structure differs between colelablinux and nm3)
#This version also includes parms that can run portions of the fMRIVol pipeline (the most time-intensive one) if this gets aborted e.g. from pre-emption

#REQUIREMENTS:
#Transfer: input MRI data in BIDS format (e.g. mriqc_BIDS), HCP_v2_prereqs directory (contains program versions appropriate for HCP pipeline AND modified matlab functions within FIX dir)
#Matlab installation and modify PATH; **PENDING R installation with appropriate packages/versions (see FSL_FIX readme) and modify PATH, modify settings.sh in FIX dir**
#Setup and transfer: Step1 pipeline script (this one), Step1a path sourcing script (to various toolboxes), opts.shlib command line script (should not need to modify), Step1b nm3 batch script (how the function is actually executed on nm3; all these 3 need to be in same scripts directory

# Cole Lab incorporation of HCP Pipelines for use in analysis of IndivRITL data set (RUBIC protocol)
# Will quit script if there is an error detected
#set -e
#########################
source opts.shlib # command line option functions
########################


show_usage() {
    cat <<EOF

This script serves as a function to run an individual subject through the HCP Pipeline, along with other necessary commands to fit their pipeline to the IndivRITL analysis.

hcp_preprocessing_func.sh

Usage: hcp_preprocessing_func.sh [options]

    --server=<servername>       (required) The name of the server you are running this from (i.e., "colelablinux", "nm3")
    --subj=<subjectID>          (required) The subject number
    --dirSetUp=<"true">         (optional) Input "true", if you want to set up directory structure for this subject
    --anatomicalRecon=<"true">  (optional) Input "true", if you want to run a DICOM reconstruction of the anatomicals
    --epiRecon=<"true">         (optional) Input "true", if you want to run a DICOM reconstruction of the EPIs
    --preFS=<"true">            (optional) Input "true", if you want to run the PreFS HCP Node
    --FS=<"true">               (optional) Input "true", if you want to run the FS HCP Node
    --postFS=<"true">           (optional) Input "true", if you want to run the Post FS HCP Node
    --fmriVol=<"true">          (optional) Input "true", if you want to run the fMRIVolume processing of the HCP Node
    --fmriVol_prac=<"true">		(optional) Input "true", if you want to run the fMRI Volume processing *for specific run numbers of PRACTICE task*
    --fmriVol_prac_run=<"Run Number"> (optional) Input a run number integer, which will run fMRI Volume processing for that run number onwards for PRACTICE
    --fmriVol_test=<"true">		(optional) Input "true", if you want to run the fMRI Volume processing *for specific run numbers of TEST task*
    --fmriVol_test_run=<"Run Number"> (optional) Input a run number integer, which will run fMRI Volume processing for that run number onwards for TEST
    --fmriSurf=<"true">         (optional) Input "true", if you want to run the fMRISurface processing of the HCP node
    --restFix=<"true">			(optional) Input "true", if you want to run FIX-ICA on the rest data (to generate inputs for MSMAll)
    --msmAll=<"true">			(optional) Input "true", if you want to run the MSMAll multimodal surface alignment HCP node (requires MSMsulc RegName in postFS)
    --dedriftResample=<"true">	(optional) Input "true", if you wan tot run the DeDriftAndResample HCP node (where MSMAll is actually applied)
    --diff=<"true">             (optional) Input "true", if you want to run the Diffusion processing of the HCP node (NOTE: As of 1/29/15, this is still not functioning properly!)
    --createMasks=<"true">           (optional) Input "true", if you want to create whole brain, gray, white and ventricle masks for fMRI volume processed data (non-HCP node)
    --concatenateRuns=<"true">  (optional) Input "true", if you want to concatenate all task runs into analysis directory, titled Task_allruns.nii.gz
    --tsExtract=<"true">        (optional) Input "true", if you want to extract the time series for whole brain, white matter, ventricle time series of both Rest1.nii.gz and Task_allruns.nii.gz. Outputs timeseries into subject's analysis directory
EOF
    exit 1
}


opts_ShowVersionIfRequested $@

if opts_CheckForHelpRequest $@; then
    show_usage
fi

#Input Variables
server=`opts_GetOpt1 "--server" $@`
subj=`opts_GetOpt1 "--subj" $@`
dirSetUp=`opts_GetOpt1 "--dirSetUp" $@`
anatomicalRecon=`opts_GetOpt1 "--anatomicalRecon" $@`
epiRecon=`opts_GetOpt1 "--epiRecon" $@`
preFS=`opts_GetOpt1 "--preFS" $@`
FS=`opts_GetOpt1 "--FS" $@`
postFS=`opts_GetOpt1 "--postFS" $@`
fmriVol=`opts_GetOpt1 "--fmriVol" $@`
fmriVol_prac=`opts_GetOpt1 "--fmriVol_prac" $@`
fmriVol_prac_run=`opts_GetOpt1 "--fmriVol_prac_run" $@`
fmriVol_test=`opts_GetOpt1 "--fmriVol_test" $@`
fmriVol_test_run=`opts_GetOpt1 "--fmriVol_test_run" $@`
fmriSurf=`opts_GetOpt1 "--fmriSurf" $@`
restFix=`opts_GetOpt1 "--restFix" $@`
msmAll=`opts_GetOpt1 "--msmAll" $@`
dedriftResample=`opts_GetOpt1 "--dedriftResample" $@`
diff=`opts_GetOpt1 "--diff" $@`
createMasks=`opts_GetOpt1 "--createMasks" $@`
concatenateRuns=`opts_GetOpt1 "--concatenateRuns" $@`
tsExtract=`opts_GetOpt1 "--tsExtract" $@`

#****RM NOTE - modify the 'amarel' section in the if statement  below to suit input/output dirs for your project; make sure that the analysis parms in the sections that follow
#are appropriate for your sequence (consult HCP github, or the example scripts for each module in HCP_v2_prereqs/HCP_Pipelines_v3_25_1/Examples/Scripts if unsure); make sure that
#parts of the volume/surface/fix/MSMall/dedriftresample pipelines that generate filenames for each subject contain appropriate *strings* for how you named your sequences

#########################
# Set up HCP Environment - This is the customized environment for ColeLabMac Server
# Figure out which server is being used (mac or linux)
# *modified colelablinux loop: HCPPipe=/usr/local/HCP_Pipelines_v2; EnvScript=${HCPPipe}/Examples/Scripts/SetUpHCPPipeline.sh
if [ -z "$server" ]; then
    echo "Missing required option. Indicate which server you're using! Exiting script..."
    exit
elif [ "$server" == "colelabmac" ]; then
    HCPPipe=/Applications/HCP_Pipelines/Pipelines-3.4.1
    EnvScript=${HCPPipe}/Examples/Scripts/SetUpHCPPipeline_Custom.sh
elif [ "$server" == "colelablinux" ]; then
    HCPPipe=/projects3/CPRO2_learning/HCP_v2_prereqs/HCP_Pipelines_v3_25_1
    EnvScript=/projects3/CPRO2_learning/docs/scripts/Step1a_SetUpHCPPipeline_msmall.sh
elif [ "$server" == "nm3" ]; then
    HCPPipe=/home/rdm146/HCP_v2_prereqs/HCP_Pipelines_v3_25_1
    EnvScript=/home/work/rdm146/projects/CPRO2_learning/docs/scripts/Step1a_SetUpHCPPipeline_msmall_nm3.sh
elif [ "$server" == "amarel" ]; then
   HCPPipe=/projects/f_keanebp/HCP_v2_prereqs/HCP_Pipelines_v3_25_1
  EnvScript=/projects/f_keanebp/NeuralMechanisms/docs/scripts/Step1a_SetUpHCPPipeline_msmall_amarel.sh
	basedir="/projects/f_keanebp/NeuralMechanisms/data/${subj}"
	datadir="${basedir}/preprocessed/"
	if [ ! -e $datadir ]; then mkdir -p $datadir; fi
	subjdir=${datadir}
	if [ ! -e $subjdir ]; then mkdir -p $subjdir; fi
	unprocesseddir="${basedir}/"
	if [ ! -e $unprocesseddir ]; then mkdir $unprocesseddir; fi
	analysisdir=${subjdir}/analysis
	if [ ! -e $analysisdir ]; then mkdir $analysisdir; fi
	subjmaskdir=${subjdir}/masks
	if [ ! -e $subjmaskdir ]; then mkdir $subjmaskdir; fi
elif [ "$server" == "amarel_mnt" ]; then
   # Updated by bpk on 7/30/2018
    HCPPipe=/home/rdm146/HCP_v2_prereqs/HCP_Pipelines_v3_25_1
    EnvScript=/home/rdm146/projects/CPRO2_learning/docs/scripts/Step1a_SetUpHCPPipeline_msmall_amarel.sh
	basedir="/mnt/scratch/rdm146/CPRO2_learning_${subj}"
        datadir="${basedir}/output/"
        if [ ! -e $datadir ]; then mkdir -p $datadir; fi
        subjdir=${datadir}
        if [ ! -e $subjdir ]; then mkdir -p $subjdir; fi
        unprocesseddir="${basedir}/input/"
        if [ ! -e $unprocesseddir ]; then mkdir $unprocesseddir; fi
        analysisdir=${subjdir}/analysis
        if [ ! -e $analysisdir ]; then mkdir $analysisdir; fi
        subjmaskdir=${subjdir}/masks
        if [ ! -e $subjmaskdir ]; then mkdir $subjmaskdir; fi
fi
# Set up HCP Pipeline Environment
. ${EnvScript}
########################

#########################
# HCP Conventions and Parameters - shouldn't need to edit this
#PostFreeSurfer input file names
#Input Variables
SurfaceAtlasDIR="${HCPPIPEDIR_Templates}/standard_mesh_atlases"
GrayordinatesSpaceDIR="${HCPPIPEDIR_Templates}/91282_Greyordinates"
GrayordinatesResolutions="2" #Usually 2mm, if multiple delimit with @, must already exist in templates dir
HighResMesh="164" #Usually 164k vertices
LowResMeshes="32" #Usually 32k vertices, if multiple delimit with @, must already exist in templates dir
SubcorticalGrayLabels="${HCPPIPEDIR_Config}/FreeSurferSubcorticalLabelTableLut.txt"
FreeSurferLabels="${HCPPIPEDIR_Config}/FreeSurferAllLut.txt"
ReferenceMyelinMaps="${HCPPIPEDIR_Templates}/standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii"
# RegName="MSMSulc" #MSMSulc is recommended, if binary is not available use FS (FreeSurfer)
RegName="MSMSulc"
# set bandpass filter used in FIX-ICA in seconds (2000 corresponds to detrending i.e. very lenient highpass 1/2000Hz)
bandpass=2000
#########################

#########################
# Data and scan parameters
procAllFuncScans="y" # Process all functional scans
# Notes, bpk: Everything after and including S31, C23, and B27, use 75 rather than 69 for the echo spacing values (DwellTime)
listOfSubjectsNoViscomp="sub-C10"
listOfSubjectsNoRetino=""
listOfSubjectsNoContour=""
listOfSubjectsNoEbb=""
listOfSubjectsOldProtocol="sub-B04 sub-B06 sub-B10 sub-B14 sub-B21 sub-B23 sub-C02 sub-C05 sub-C10 sub-C13 sub-C14 sub-C15 sub-C22 sub-S06 sub-S10 sub-S12 sub-S16 sub-S20 sub-S25"
SmoothingFWHM="2"
if [[ $listOfSubjectsOldProtocol == *"${subj}"* ]]; then
    DwellTime_SE="0.00069" # the dwell time or echo spacing of the SE FieldMaps (see protocol)
    DwellTime_fMRI="0.00069" # the dwell time or echo spacing of the fMRI multiband sequence (see protocol)
    echo "Subject ${subj} uses the old protocol."
elif [[ $listOfSubjectsOldProtocol != *"${subj}"* ]]; then
    DwellTime_SE="0.00075" # the dwell time or echo spacing of the SE FieldMaps (see protocol)
    DwellTime_fMRI="0.00075" # the dwell time or echo spacing of the fMRI multiband sequence (see protocol)
    echo "Subject ${subj} uses the new protocol."
fi
T1wSampleSpacing="0.0000074" # This parameter can be found at DICOM field (0019,1018) (use command `dicom_hdr *001.dcm | grep "0019 1018"`
T2wSampleSpacing="0.0000021" # This parameter can be found at DICOM field (0019,1018) (use command `dicom_hdr *001.dcm | grep "0019 1018"`
# Default parameters (should not need to be changed for this study
fmrires="2.4"
brainsize="150"
seunwarpdir="y" # AP/PA is Y
unwarpdir="-y" # unwarp direction for A >> P phase encoded data is -y (think about how anterior to posterior coordinate direction is -y). It follows that unwarpdir for P >> A collected data would be "y", but we have not collected this data for the IndivRITL study.
#numTRsPerTaskRun=281 # updated by bpk for VisComp; apparently this is not needed below.

# Anatomical templates for this study (MNI templates)
t1template="${HCPPipe}/global/templates/MNI152_T1_0.8mm.nii.gz"
t1template2mm="${HCPPipe}/global/templates/MNI152_T1_2mm.nii.gz"
t1templatebrain="${HCPPipe}/global/templates/MNI152_T1_0.8mm_brain.nii.gz"
t2template="${HCPPipe}/global/templates/MNI152_T2_0.8mm.nii.gz"
t2templatebrain="${HCPPipe}/global/templates/MNI152_T2_0.8mm_brain.nii.gz"
t2template2mm="${HCPPipe}/global/templates/MNI152_T2_2mm.nii.gz"
templatemask="${HCPPipe}/global/templates/MNI152_T1_0.8mm_brain_mask.nii.gz"
template2mmmask="${HCPPipe}/global/templates/MNI152_T1_2mm_brain_mask_dil.nii.gz"


#########################
# Start first node of HCP Pipeline: PreFreeSurferPipeline
#fmaps needs to be from run2, cos that's when anatomicals were collected
#**note option to also pass gre fieldmaps using --fmapmag and --fmapphase
if [ -z "$preFS" ]; then
    echo "Skipping PreFS HCP Node"
elif [ $preFS = true ]; then
     echo ${HCPPipe}/PreFreeSurfer/PreFreeSurferPipeline.sh
    # changed rec-1 to rec-2 and changed SEPhase variables to run-1, bpk, 07/2018
    ${HCPPipe}/PreFreeSurfer/PreFreeSurferPipeline.sh \
        --path="${datadir}" \
        --subject="${subj}" \
        --t1="${unprocesseddir}/anat/${subj}_rec-1_T1w.nii.gz" \
        --t2="${unprocesseddir}/anat/${subj}_rec-1_T2w.nii.gz" \
        --t1template="${t1template}" \
        --t1templatebrain="${t1templatebrain}" \
        --t1template2mm="${t1template2mm}" \
        --t2template="${t2template}" \
        --t2templatebrain="$t2templatebrain" \
        --t2template2mm="$t2template2mm" \
        --templatemask="$templatemask" \
        --template2mmmask="$template2mmmask" \
        --brainsize="${brainsize}" \
        --fmapmag="NONE" \
        --fnirtconfig="${HCPPipe}/global/config/T1_2_MNI152_2mm.cnf" \
        --SEPhaseNeg="${unprocesseddir}/fmap/${subj}_dir-AP_run-1_epi.nii.gz" \
        --SEPhasePos="${unprocesseddir}/fmap/${subj}_dir-PA_run-1_epi.nii.gz" \
        --echospacing="$DwellTime_SE" \
        --seunwarpdir="${seunwarpdir}" \
        --t1samplespacing="$T1wSampleSpacing" \
        --t2samplespacing="$T2wSampleSpacing" \
        --unwarpdir="z" \
        --grdcoeffs="NONE" \
        --avgrdcmethod="TOPUP" \
        --topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf" \
        --printcom=""
fi

#########################
# Start second node of HCP Pipeline: FreeSurferPipeline
if [ -z "$FS" ]; then
    echo "Skipping Freesurfer HCP node"
elif [ $FS = true ]; then
    # limit number of threads of FS
    export OMP_NUM_THREADS=3 #
    ${HCPPipe}/FreeSurfer/FreeSurferPipeline.sh \
        --subject="${subj}" \
        --subjectDIR="${datadir}/${subj}/T1w" \
        --t1="${datadir}/${subj}/T1w/T1w_acpc_dc_restore.nii.gz" \
        --t1brain="${datadir}/${subj}/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
        --t2="${datadir}/${subj}/T1w/T2w_acpc_dc_restore.nii.gz"
fi

#########################
# Start third node of HCP Pipeline: PostFreeSurferPipeline
if [ -z "$postFS" ]; then
    echo "Skipping PostFS HCP node"
elif [ $postFS = true ]; then

    ${HCPPipe}/PostFreeSurfer/PostFreeSurferPipeline.sh \
        --path="${datadir}" \
        --subject="${subj}" \
        --surfatlasdir="$SurfaceAtlasDIR" \
        --grayordinatesdir="$GrayordinatesSpaceDIR" \
        --grayordinatesres="$GrayordinatesResolutions" \
        --hiresmesh="$HighResMesh" \
        --lowresmesh="$LowResMeshes" \
        --subcortgraylabels="$SubcorticalGrayLabels" \
        --freesurferlabels="$FreeSurferLabels" \
        --refmyelinmaps="$ReferenceMyelinMaps" \
        --regname="$RegName" \
        --printcom=""
fi

#########################
# Start fourth node of HCP Pipeline: GenericfMRIVolumeProcessing
#*one block for each of the four tasks: Pracrest, Practask, Testrest, Testtask
if [ -z "$fmriVol" ]; then
    echo "Skipping fMRIVolumeProcessing node"
elif [ $fmriVol = true ]; then
    # need to iterate through each rest scan, and then each task scan
    # note in this version you have to specify --biascorrection method; following parm is taken from the HCP volume example script; previous lab version of the pipeline might have used LEGACY, but SEBASED is better
    BiasCorrection="SEBASED" #NONE, LEGACY, or SEBASED: LEGACY uses the T1w bias field, SEBASED calculates bias field from spin echo images (which requires TOPUP distortion correction)

    if [[ $procAllFuncScans == "y" ]]; then
      ## Rest scan, updated by bpk on 7/31/18, using run 2 of field maps for this part only.
    	echo "Running fMRI Volume processing on Rest scan"
    	fmriname="Task_Rest"
    	fmritcs="${unprocesseddir}/func/${subj}_task-rest_run-1_bold.nii.gz"
    	fmriscout="${unprocesseddir}/func/${subj}_task-rest_run-1_sbref.nii.gz"
    	fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-2_epi.nii.gz"
    	fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-2_epi.nii.gz"
    	${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
    		--path="${datadir}" \
    		--subject="${subj}" \
    		--fmriname="${fmriname}" \
    		--fmritcs="${fmritcs}" \
    		--fmriscout="${fmriscout}" \
    		--SEPhaseNeg="${fmap_neg_ap}" \
    		--SEPhasePos="${fmap_pos_pa}" \
    		--fmapmag="NONE" \
    		--fmapphase="NONE" \
    		--echospacing="$DwellTime_fMRI" \
    		--echodiff="NONE" \
    		--unwarpdir="${unwarpdir}" \
    		--fmrires="$fmrires" \
    		--dcmethod="TOPUP" \
    		--gdcoeffs="NONE" \
    		--printcom="" \
    		--biascorrection=${BiasCorrection} \
    		--topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"

        ## Viscomp Task Practice Scans, updated by bpk, 7/30/201
        pushd ${unprocesseddir}/func/
        numPracRuns=`ls *viscomp*bold.nii.gz | wc -l` # spits out 4.
        echo "Viscomp task runs for $subj are $numPracRuns"
        popd
        for ((i=1;i<=${numPracRuns};++i)); do
            # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
            echo "Running fMRI Volume processing on Viscomp scan #${i}"
            fmriname="Task_Viscomp${i}" # Changed by bpk, 7/31/2018
            fmritcs="${unprocesseddir}/func/${subj}_task-viscomp_run-${i}_bold.nii.gz"
    		    fmriscout="${unprocesseddir}/func/${subj}_task-viscomp_run-${i}_sbref.nii.gz"
    		    fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-1_epi.nii.gz"
    		    fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-1_epi.nii.gz"
            ${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
                --path="${datadir}" \
                --subject="${subj}" \
                --fmriname="${fmriname}" \
                --fmritcs="${fmritcs}" \
                --fmriscout="${fmriscout}" \
                --SEPhaseNeg="${fmap_neg_ap}" \
                --SEPhasePos="${fmap_pos_pa}" \
                --fmapmag="NONE" \
                --fmapphase="NONE" \
                --echospacing="$DwellTime_fMRI" \
                --echodiff="NONE" \
                --unwarpdir="${unwarpdir}" \
                --fmrires="$fmrires" \
                --dcmethod="TOPUP" \
                --gdcoeffs="NONE" \
                --printcom="" \
                --biascorrection=${BiasCorrection} \
                --topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"
        done

        ## Retinotopy Scans
        #first need to generate number of scans (given aborted runs etc)
        pushd ${unprocesseddir}/func/
        numTestRuns=`ls *retino*bold.nii.gz | wc -l`
        echo "Retino task runs for $subj are $numTestRuns"
        popd
        for ((i=1;i<=${numTestRuns};++i)); do
            # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
            echo "Running fMRI Volume processing on Retino Task scan #${i}"
            fmriname="Task_Retino${i}"
            fmritcs="${unprocesseddir}/func/${subj}_task-retino_run-${i}_bold.nii.gz"
    		    fmriscout="${unprocesseddir}/func/${subj}_task-retino_run-${i}_sbref.nii.gz"
    		    fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-1_epi.nii.gz"
    		    fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-1_epi.nii.gz"
            ${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
                --path="${datadir}" \
                --subject="${subj}" \
                --fmriname="${fmriname}" \
                --fmritcs="${fmritcs}" \
                --fmriscout="${fmriscout}" \
                --SEPhaseNeg="${fmap_neg_ap}" \
                --SEPhasePos="${fmap_pos_pa}" \
                --fmapmag="NONE" \
                --fmapphase="NONE" \
                --echospacing="$DwellTime_fMRI" \
                --echodiff="NONE" \
                --unwarpdir="${unwarpdir}" \
                --fmrires="$fmrires" \
                --dcmethod="TOPUP" \
                --gdcoeffs="NONE" \
                --printcom="" \
                --biascorrection=${BiasCorrection} \
                --topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"
        done
    fi # done with if-statement, processing all functional scans

    ## Contour Scans
    #first need to generate number of scans (given aborted runs etc)
    pushd ${unprocesseddir}/func/
    numTestRuns=`ls *contour*bold.nii.gz | wc -l`
    echo "Contour task runs for $subj are $numTestRuns"
    popd
    for ((i=1;i<=${numTestRuns};++i)); do
        # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
        echo "Running fMRI Volume processing on Contour Task scan #${i}"
        fmriname="Task_Contour${i}"
        fmritcs="${unprocesseddir}/func/${subj}_task-contour_run-${i}_bold.nii.gz"
		    fmriscout="${unprocesseddir}/func/${subj}_task-contour_run-${i}_sbref.nii.gz"
		    fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-2_epi.nii.gz"
		    fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-2_epi.nii.gz"
        ${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
            --path="${datadir}" \
            --subject="${subj}" \
            --fmriname="${fmriname}" \
            --fmritcs="${fmritcs}" \
            --fmriscout="${fmriscout}" \
            --SEPhaseNeg="${fmap_neg_ap}" \
            --SEPhasePos="${fmap_pos_pa}" \
            --fmapmag="NONE" \
            --fmapphase="NONE" \
            --echospacing="$DwellTime_fMRI" \
            --echodiff="NONE" \
            --unwarpdir="${unwarpdir}" \
            --fmrires="$fmrires" \
            --dcmethod="TOPUP" \
            --gdcoeffs="NONE" \
            --printcom="" \
            --biascorrection=${BiasCorrection} \
            --topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"
    fi

    ## EbbLoc Scan
    ## Rest scan, updated by bpk on 7/31/18, using run 2 of field maps for this part only.
  	echo "Running fMRI Volume processing on EbbLoc scan"
  	fmriname="Task_EbbLoc"
  	fmritcs="${unprocesseddir}/func/${subj}_task-ebbloc_run-1_bold.nii.gz"
  	fmriscout="${unprocesseddir}/func/${subj}_task-ebbloc_run-1_sbref.nii.gz"
  	fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-2_epi.nii.gz"
  	fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-2_epi.nii.gz"
  	${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
  		--path="${datadir}" \
  		--subject="${subj}" \
  		--fmriname="${fmriname}" \
  		--fmritcs="${fmritcs}" \
  		--fmriscout="${fmriscout}" \
  		--SEPhaseNeg="${fmap_neg_ap}" \
  		--SEPhasePos="${fmap_pos_pa}" \
  		--fmapmag="NONE" \
  		--fmapphase="NONE" \
  		--echospacing="$DwellTime_fMRI" \
  		--echodiff="NONE" \
  		--unwarpdir="${unwarpdir}" \
  		--fmrires="$fmrires" \
  		--dcmethod="TOPUP" \
  		--gdcoeffs="NONE" \
  		--printcom="" \
  		--biascorrection=${BiasCorrection} \
  		--topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"

      ## EbbReg Scans
      #first need to generate number of scans (given aborted runs etc)
      pushd ${unprocesseddir}/func/
      numTestRuns=`ls *ebbreg*bold.nii.gz | wc -l`
      echo "EbbReg task runs for $subj are $numTestRuns"
      popd
      for ((i=1;i<=${numTestRuns};++i)); do
          # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
          echo "Running fMRI Volume processing on EbbReg Task scan #${i}"
          fmriname="Task_EbbReg${i}"
          fmritcs="${unprocesseddir}/func/${subj}_task-ebbreg_run-${i}_bold.nii.gz"
  		    fmriscout="${unprocesseddir}/func/${subj}_task-ebbreg_run-${i}_sbref.nii.gz"
  		    fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-2_epi.nii.gz"
  		    fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-2_epi.nii.gz"
          ${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
              --path="${datadir}" \
              --subject="${subj}" \
              --fmriname="${fmriname}" \
              --fmritcs="${fmritcs}" \
              --fmriscout="${fmriscout}" \
              --SEPhaseNeg="${fmap_neg_ap}" \
              --SEPhasePos="${fmap_pos_pa}" \
              --fmapmag="NONE" \
              --fmapphase="NONE" \
              --echospacing="$DwellTime_fMRI" \
              --echodiff="NONE" \
              --unwarpdir="${unwarpdir}" \
              --fmrires="$fmrires" \
              --dcmethod="TOPUP" \
              --gdcoeffs="NONE" \
              --printcom="" \
              --biascorrection=${BiasCorrection} \
              --topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"
    done
fi


#########################
# Start fifth node of HCP Pipeline: GenericfMRISurfaceProcessing
#*one block for each of the four tasks: Pracrest, Practask, Testrest, Testtask
if [ -z "$fmriSurf" ]; then
    echo "Skipping fMRI Surface Processing node"
elif [ $fmriSurf = true ]; then

    if [[ $procAllFuncScans == "y" ]]; then
        # set input name
        fmriname="Task_Rest"
        ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
       		--path="${datadir}" \
            --subject="${subj}" \
            --fmriname="${fmriname}" \
            --fmrires="$fmrires" \
            --lowresmesh="$LowResMeshes" \
            --smoothingFWHM=$SmoothingFWHM \
            --grayordinatesres="$GrayordinatesResolutions" \
            --regname=$RegName

        ## Task Viscomp
        #first need to generate number of scans (given aborted runs etc)
        pushd ${unprocesseddir}/func/
        numPracRuns=`ls *viscomp*bold.nii.gz | wc -l`
        echo "Viscomp Task runs for $subj are $numPracRuns"
        popd
        for ((i=1;i<=${numPracRuns};++i)); do
            # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
            # set input name
            fmriname="Task_Viscomp${i}"
            ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
                --path="${datadir}" \
                --subject="${subj}" \
                --fmriname="${fmriname}" \
                --fmrires="$fmrires" \
                --lowresmesh="$LowResMeshes" \
                --smoothingFWHM=$SmoothingFWHM \
                --grayordinatesres="$GrayordinatesResolutions" \
                --regname=$RegName
        done

        ## Task Retino
        #first need to generate number of scans (given aborted runs etc)
        pushd ${unprocesseddir}/func/
        numTestRuns=`ls *retino*bold.nii.gz | wc -l`
        echo "Retino Task runs for $subj are $numTestRuns"
        popd
        for ((i=1;i<=${numTestRuns};++i)); do
            # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
            # set input name
            fmriname="Task_Retino${i}"
            ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
                --path="${datadir}" \
                --subject="${subj}" \
                --fmriname="${fmriname}" \
                --fmrires="$fmrires" \
                --lowresmesh="$LowResMeshes" \
                --smoothingFWHM=$SmoothingFWHM \
                --grayordinatesres="$GrayordinatesResolutions" \
                --regname=$RegName
        done
    fi # end of if-statement, processing all functional scans.


    ## Task Contour
    #first need to generate number of scans (given aborted runs etc)
    pushd ${unprocesseddir}/func/
    numTestRuns=`ls *contour*bold.nii.gz | wc -l`
    echo "Contour Task runs for $subj are $numTestRuns"
    popd
    for ((i=1;i<=${numTestRuns};++i)); do
        # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename sub-.
        # set input name
        fmriname="Task_Contour${i}"
        ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
            --path="${datadir}" \
            --subject="${subj}" \
            --fmriname="${fmriname}" \
            --fmrires="$fmrires" \
            --lowresmesh="$LowResMeshes" \
            --smoothingFWHM=$SmoothingFWHM \
            --grayordinatesres="$GrayordinatesResolutions" \
            --regname=$RegName
    done

    # set input name
    fmriname="Task_EbbLoc"
    ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
   		  --path="${datadir}" \
        --subject="${subj}" \
        --fmriname="${fmriname}" \
        --fmrires="$fmrires" \
        --lowresmesh="$LowResMeshes" \
        --smoothingFWHM=$SmoothingFWHM \
        --grayordinatesres="$GrayordinatesResolutions" \
        --regname=$RegName

      ## Task EbbReg
      #first need to generate number of scans (given aborted runs etc)
      pushd ${unprocesseddir}/func/
      numTestRuns=`ls *ebbreg*bold.nii.gz | wc -l`
      echo "EbbReg Task runs for $subj are $numTestRuns"
      popd
      for ((i=1;i<=${numTestRuns};++i)); do
          # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename sub-.
          # set input name
          fmriname="Task_EbbReg${i}"
          ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
              --path="${datadir}" \
              --subject="${subj}" \
              --fmriname="${fmriname}" \
              --fmrires="$fmrires" \
              --lowresmesh="$LowResMeshes" \
              --smoothingFWHM=$SmoothingFWHM \
              --grayordinatesres="$GrayordinatesResolutions" \
              --regname=$RegName
      done

fi

#########################
#Fix-ICA: applies spatial ICA (via MELODIC, on volume space data) and identifies artifacts via ML classifier (via FIX)
#apply *separately* to practice_rest and test_rest (given that these were acquired in separate sessions, would combine otherwise; they are sufficiently long to provide good solutions anyway, based on the Salimi-Khorshid 2014 paper)
#not using multi-run as there is only one rest run per session, and this is of sufficient length (>10min, which was used in the Salimi FIX 2014 paper)
#outputs from Fix-ICA will be used for MSMall i.e. to provide the RSN spatial maps (between area FC) that constitute 1 of the 4 modalities

#*PREREQS for running this: 1) set up Step1a_SetUP script (Fix vars), 2) modify settings.sh in the FIX dir, 3) make sure all matlab and R are installed correctly (check R package versions required in README file in FSL_DIR)
#NOTE point to hcp_fix in FSL_FIXDIR (set in Step1a) not in $HCPPipe (latter was modified to run on HCP servers, and generates errors)
if [ -z "$restFix" ]; then
    echo "Skipping fMRI rest FIX ICA module"
elif [ $restFix = true ]; then

    ## Rest Practice
    # set input name
    fmriname="Task_Rest" #Not sure if this is needed when the data will be used for MSM-all, bpk, 8/2/18
    fmri_input=${subjdir}/${subj}/MNINonLinear/Results/${fmriname}/${fmriname}.nii.gz
    ${FSL_FIXDIR}/hcp_fix ${fmri_input} ${bandpass}
fi

#########################
#MSMall: performs multi-modal surface alignment
#MSMall computes the transform, whereas DedriftAndResample corrects i.e. dedrifts the transform (not doing this currently) AND actually applies it in one step, and resamples from native to standard grayordinates

#fMRINames should only contain rest scans, as this is the module where RSN maps and visuotopic FC maps (2 of the 4 modalities used for alignment) are generated
#Details from Glasser 2016 supplementary methods (2.3): MSMall performs Weighted spatial multiple regression (over 2 rounds) is used to derive 'refined subject-specific spatial maps', f
#or RSN network maps; mapping group ICA template maps (representing likely canonical RSNs) with subject ICA FIX maps in a refined way.

#After performing the weighted regression to generate subject-specific RSN info, MSMall was performed as follows (Glasser 2016 supplement, 2.4):
#The following 44 maps were jointly used in the final iteration of MSMAll: 34 RSNs (from Fix ICA), the subject's myelin map, eight V1-based rfMRI visotopic regressor maps (see below #4.4; based on FC gradients), and a binary
#non-cortical medial wall ROI (i.e. the region of the surface mesh that does not contain neocortical grey matter). Group versions of all these maps served as the multimodal registration target.
#TAKE HOME MESSAGE: Alignment of these subject-specific maps with their corresponding group templates (all represented as spherical surfaces) was performed using modified versions of the MSM algorithm reported in Robinson et al 2014
if [ -z "$msmAll" ]; then
    echo "Skipping MSMAll multi-modal surface alignment"
elif [ $msmAll = true ]; then
	#@ separated list of scans to perform MSMall alignment on - performs concatenation of rest scans (across runs/sessions), and generation of subject RSN ICA (after regression with group RSN templates)
	#fMRINames="Rest_Practice@Rest_Test"
	fMRINames="Task_Rest"
	#name to give to concatenated single subject "scan">
	OutfMRIName="Task_Rest_All"
	#Identification for FIX cleaned dtseries to use
	#The dense timeseries files used will be named <fmri_name>_<fmri_proc_string>.dtseries.nii where
    #fmri_name> is each of the fMRIs specified in the <fMRI Names> list and <fmri_proc_string> is this specified value
	fMRIProcSTRING="_Atlas_hp${bandpass}_clean"
	#path to group templates generated in Glasser 2016 from group ICA - weighted regression is used to generate subject-specific maps, which are then used for MSMall alignment (in 2 steps)
	MSMAllTemplates="${HCPPIPEDIR}/global/templates/MSMAll"
	#name to give output registration
	RegName_MSMall="MSMAll_InitalReg"
	#HighResMesh=${HighResMesh}
	#LowResMesh=${LowResMeshes}
	#input registration - from PostFreesurfer step (i.e. MSMSulc)
	InRegName=${RegName}
	MatlabMode="1" #Mode=0 compiled Matlab, Mode=1 interpreted Matlab; *can run both but going with interpreted for now (compiled allows Fix to run without a matlab license)
	${HCPPipe}/MSMAll/MSMAllPipeline.sh \
	  --path=${datadir} \
	  --subject=${subj} \
	  --fmri-names-list=${fMRINames} \
	  --output-fmri-name=${OutfMRIName} \
	  --high-pass=${bandpass} \
	  --fmri-proc-string=${fMRIProcSTRING} \
	  --msm-all-templates=${MSMAllTemplates} \
	  --output-registration-name=${RegName_MSMall} \
	  --high-res-mesh=${HighResMesh} \
	  --low-res-mesh=${LowResMeshes} \
	  --input-registration-name=${InRegName} \
	  --matlab-run-mode=${MatlabMode}
fi


#########################
#Dedrift and Resample: computes group 'drift' in MSM alignment and corrects for this (optional step, NOT applying for now as prereq MSMRemoveGroupDrift code is hard to parse)
#then actually applies the MSMall transform + resamples to specified func files, and specifies whether FIX-ICA is applied to resulting surface files
#check that it is feasible to provide same filenames to rfmrinames and tfmrinames - enables comparison of FIX-ICA for rest and task - NOW not running ICA FIX on rest or task, so only fill out tfmrinames
#check that RegName_out is correct below (based on previous MSMall output)

if [ -z "$dedriftResample" ]; then
    echo "Skipping DeDrift and Resample (apply MSMAll)"
elif [ "$dedriftResample" == "true" ]; then
	#HighResMesh="164"
	#LowResMesh="32"

	#String corresponding to the MSMAll or other registration sphere name e.g. ${Subject}.${Hemisphere}.sphere.${RegName}.native.surf.gii
	RegName_out="MSMAll_InitalReg_2_d40_WRN"
	#Path to the spheres output by the MSMRemoveGroupDrift pipeline or NONE
	#from example"${HCPPIPEDIR}/global/templates/MSMAll/DeDriftingGroup.L.sphere.DeDriftMSMAll.164k_fs_LR.surf.gii@${HCPPIPEDIR}/global/templates/MSMAll/DeDriftingGroup.R.sphere.DeDriftMSMAll.164k_fs_LR.surf.gii"
	DeDriftRegFiles="NONE"
	#String corresponding to the output name of the concatenated registration i.e. the dedrifted registration
	#*given that we are not running the group dedrift (Glasser says this is optional and the code is not easy to parse), seems like this should equal RegName_out (based on looking at the DeDrift shell script)
	ConcatRegName=${RegName_out}
	#delimited map name strings corresponding to maps that are not myelin maps e.g. sulc curvature corrThickness thickness
	Maps="sulc@curvature@corrThickness@thickness"
	#delimited map name strings corresponding to myelin maps e.g. MyelinMap SmoothedMyelinMap) No _BC, this will be reapplied
	MyelinMaps="MyelinMap@SmoothedMyelinMap"
  if [[ $procAllFuncScans == "y" ]]; then
    	#generate run names for Task_Retino and Task_Viscomp
        pushd ${unprocesseddir}/func/
        numViscompRuns=`ls *viscomp*bold.nii.gz | wc -l`
        echo " Task Viscomp runs for $subj are $numViscompRuns"
        popd

        pushd ${unprocesseddir}/func/
        numRetinoRuns=`ls *retino*bold.nii.gz | wc -l`
        echo "Test Retino task runs for $subj are $numRetinoRuns"
        popd
  fi

  pushd ${unprocesseddir}/func/
  numContourRuns=`ls *contour*bold.nii.gz | wc -l`
  echo "Test Contour task runs for $subj are $numContourRuns"
  popd

  pushd ${unprocesseddir}/func/
  numEbbLocRuns=`ls *ebbloc*bold.nii.gz | wc -l`
  echo "Test EbbLoc task runs for $subj are $numEbbLocRuns"
  popd

  pushd ${unprocesseddir}/func/
  numEbbRegRuns=`ls *ebbreg*bold.nii.gz | wc -l`
  echo "Test EbbReg task runs for $subj are $numEbbRegRuns"
  popd

  if [[ $procAllFuncScans == "y" ]]; then
    	run_names="Task_Viscomp1" #initial string
        for ((i=2;i<=${numViscompRuns};++i)); do
    		run_names="${run_names}@Task_Viscomp${i}"
        done

    	run_names="${run_names}@Task_Retino1"
    	  for ((i=2;i<=${numRetinoRuns};++i)); do
    		run_names="${run_names}@Task_Retino${i}"
        done
  fi

  if [[ $procAllFuncScans == "y" ]]; then
      run_names="${run_names}@Task_Contour1" # this statement is correct!
  elif [[ $procAllFuncScans == "n" ]]; then
    run_names="Task_Contour1"
  fi
  for ((i=2;i<=${numContourRuns};++i)); do
    run_names="${run_names}@Task_Contour${i}"
  done

  run_names="${run_names}@Task_EbbLoc"
#    for ((i=2;i<=${numRetinoRuns};++i)); do
#    run_names="${run_names}@Task_EbbLoc${i}"
#    done

  run_names="${run_names}@Task_EbbReg1"
    for ((i=2;i<=${numEbbRegRuns};++i)); do
    run_names="${run_names}@Task_EbbReg${i}"
    done

	#*@ delimited fMRIName strings corresponding to maps that will have ICA+FIX reapplied to them (could be either rfMRI or tfMRI). If none are to be used, specify "NONE".
	#rfMRINames="rfMRI_REST1_LR rfMRI_REST1_RL rfMRI_REST2_LR rfMRI_REST2_RL" #Space delimited list or NONE
	#rfMRINames="Rest_Practice@Rest_Test@Task_Practice@Task_Test"
	rfMRINames="NONE"
	#*@ delimited fMRIName strings corresponding to maps that will not have ICA+FIX reapplied to them (likely not to be used in the future as ICA+FIX will be recommended for all fMRI data) If none are to be used, specify "NONE".
	#tfMRINames="tfMRI_WM_LR tfMRI_WM_RL tfMRI_GAMBLING_LR tfMRI_GAMBLING_RL tfMRI_MOTOR_LR tfMRI_MOTOR_RL tfMRI_LANGUAGE_LR tfMRI_LANGUAGE_RL tfMRI_SOCIAL_LR tfMRI_SOCIAL_RL tfMRI_RELATIONAL_LR tfMRI_RELATIONAL_RL tfMRI_EMOTION_LR tfMRI_EMOTION_RL" #Space delimited list or NONE
	#*not sure if run numbers need to be added for task scans?
  if [[ $procAllFuncScans == "y" ]]; then
    tfMRINames="Task_Rest@${run_names}"
  elif [[ $procAllFuncScans == "n" ]]; then
    tfMRINames="${run_names}"
  fi
	echo "tfMRINames for this subject are $tfMRINames"

	#Should equal previous grayordiantes smoothing in fMRISurface (because we are resampling from unsmoothed native mesh timeseries - already set
	MatlabMode="1" #Mode=0 compiled Matlab, Mode=1 interpreted Matlab

	${HCPPipe}/DeDriftAndResample/DeDriftAndResamplePipeline.sh \
	  --path=${datadir} \
	  --subject=${subj} \
	  --high-res-mesh=${HighResMesh} \
	  --low-res-meshes=${LowResMeshes} \
	  --registration-name=${RegName_out} \
	  --dedrift-reg-files=${DeDriftRegFiles} \
	  --concat-reg-name=${ConcatRegName} \
	  --maps=${Maps} \
	  --myelin-maps=${MyelinMaps} \
	  --rfmri-names=${rfMRINames} \
	  --tfmri-names=${tfMRINames} \
	  --smoothing-fwhm=${SmoothingFWHM} \
	  --highpass=${bandpass} \
	  --matlab-run-mode=${MatlabMode}
fi
