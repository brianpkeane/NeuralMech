#!/bin/bash
# Taku Ito
# 10/29/14

#modified by Ravi Mill for CPRO2_learning
#assumes that dicom2nifti conversion has been carried out (e.g. as part of BIDS conversion)
# 03/16/18

# Modified again by Brian Keane to preprocess data from "Neural Mechanisms" study.
# 07-08/2018
# Changed anatomical scans to reconstruction 2 (rec-2)
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
   # HCPPipe=/home/rdm146/HCP_v2_prereqs/HCP_Pipelines_v3_25_1
    HCPPipe=/projects/f_mc1689_1/HCP_v2_prereqs/HCP_Pipelines_v3_25_1 
    EnvScript=/home/keanebp/projects/NeuralMech/docs/scripts/Step1a_SetUpHCPPipeline_msmall_amarel.sh
	basedir="/scratch/keanebp/NeuralMech_${subj}"
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
## Set up Subject Directory Parameters - **Now assigned above but leaving the following commands for reference
##**make sure nm3 basedir matches with scratch dir paths in Step1b batch script
#if [ "$server" == "nm3" ]; then
#    basedir="/scratch/rdm146/CPRO2_learning_${subj}"
#elif [ "$server" == "amarel" ]; then
#    basedir="/scratch/rdm146/CPRO2_learning_${subj}"
#else
#    basedir="/projects3/CPRO2_learning/"
#fi

##root output directory
#datadir="${basedir}/output/"
#if [ ! -e $datadir ]; then mkdir -p $datadir; fi 
##subject output directory - *same as datadir for nm3 version (to enable deleting of entire subject folder, without interfering with rest of output); needed as convention for hcp modules is datadir and subjNum passed separately, except for fix module and mask modules which are set up for a combined call
#subjdir=${datadir}
#if [ ! -e $subjdir ]; then mkdir -p $subjdir; fi 

##input directory
#unprocesseddir="${basedir}/input/"
#if [ ! -e $unprocesseddir ]; then mkdir $unprocesseddir; fi 
#rawdatadir="${datadir}/rawdata/${subj}/MRI/"
#if [ ! -e $rawdatadir ]; then mkdir $rawdatadir; fi 
##subject output subdirectories - used in later GLM script
#analysisdir=${subjdir}/analysis
#if [ ! -e $analysisdir ]; then mkdir $analysisdir; fi
#subjmaskdir=${subjdir}/masks
#if [ ! -e $subjmaskdir ]; then mkdir $subjmaskdir; fi 

#########################



#########################

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
# Notes, bpk: Everything after and including S31, C23, and B27, use 75 rather than 69 for the echo spacing values (DwellTime) 
SmoothingFWHM="2.4" # Changed from 2.0, bpk, 8.6.2018
DwellTime_SE="0.00075" # the dwell time or echo spacing of the SE FieldMaps (see protocol)
DwellTime_fMRI="0.00075" # the dwell time or echo spacing of the fMRI multiband sequence (see protocol)
T1wSampleSpacing="0.0000074" # This parameter can be found at DICOM field (0019,1018) (use command `dicom_hdr *001.dcm | grep "0019 1018"`
T2wSampleSpacing="0.0000021" # This parameter can be found at DICOM field (0019,1018) (use command `dicom_hdr *001.dcm | grep "0019 1018"`
# Default parameters (should not need to be changed for this study
fmrires="2.4"
brainsize="150" 
seunwarpdir="y" # AP/PA is Y
unwarpdir="-y" # unwarp direction for A >> P phase encoded data is -y (think about how anterior to posterior coordinate direction is -y). It follows that unwarpdir for P >> A collected data would be "y", but we have not collected this data for the IndivRITL study.
numTRsPerTaskRun=281 # updated by bpk for VisComp; apparently this is not needed below.


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

#**Skipping setup of behav directories and dicom2nifti conversion
#########################
# Set up directories for behavioral data
# This is not part of the HCP pipeline
#if [ -z "$dirSetUp" ]; then
#    echo "Skipping directory set up"
#elif [ "$dirSetUp" == "true" ]; then
#    echo "Setting up directories"
#    behavdir=${rawdatadir}/behavdata
#    duncan=${behavdir}/Duncan
#    mobile=${behavdir}/MobileSurveys
#    nih=${behavdir}/NIHToolBox
#    opspan=${behavdir}/OpSpan
#    fmri=${behavdir}/fMRI_Behavioral
#    cproprac=${behavdir}/CPRO_Practice
#    if [ ! -e $behavdir ]; then mkdir $behavdir; fi 
#    if [ ! -e $duncan ]; then mkdir $duncan; fi 
#    if [ ! -e $mobile ]; then mkdir $mobile; fi 
#    if [ ! -e $nih ]; then mkdir $nih; fi 
#    if [ ! -e $opspan ]; then mkdir $opspan; fi 
#    if [ ! -e $fmri ]; then mkdir $fmri; fi 
#    if [ ! -e $cproprac ]; then mkdir $cproprac; fi
#    # Allow other users without admin privileges to edit these directories
#    chmod -R 775 ${rawdatadir}
#fi
#########################
# DICOM to NIFTI reconstruction using Freesurfer (Everything is reconstructed EXCEPT Diffusion data)
# This is NOT part of the HCP pipeline
# Reconstruct Anatomicals and Field Maps
#if [ -z "$anatomicalRecon" ]; then
    #echo "Skipping anatomical reconstruction"
#elif [ $anatomicalRecon = true ]; then


    #if [ ! -e ${unprocesseddir}/T1w ]; then mkdir ${unprocesseddir}/T1w; fi 
    #if [ ! -e ${unprocesseddir}/T2w ]; then mkdir ${unprocesseddir}/T2w; fi 
    #if [ ! -e ${unprocesseddir}/SE_Maps ]; then mkdir ${unprocesseddir}/SE_Maps; fi 
    #if [ ! -e ${unprocesseddir}/GRE_FieldMap ]; then mkdir ${unprocesseddir}/GRE_FieldMap; fi 
    
    #echo "Reconstructing T1 image"
    #mri_convert ${rawdatadir}/${T1w}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/T1w/T1w.nii.gz
    #echo "Reconstructing T2 image"
    #mri_convert ${rawdatadir}/${T2w}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/T2w/T2w.nii.gz
    #echo "Reconstructing Spin Echo Field Maps"
    #mri_convert ${rawdatadir}/${SE_AP}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/SE_Maps/SE_AP.nii.gz
    #mri_convert ${rawdatadir}/${SE_PA}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/SE_Maps/SE_PA.nii.gz
    #echo "Reconstructing Gradient Field Map"
    #mri_convert ${rawdatadir}/${GRE_FIELDMAP}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/GRE_FieldMap/GRE_FieldMap.nii.gz
    #echo "Resampling Field Map images to RPI coordinates"
    #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/SE_Maps/SE_AP.nii.gz -prefix ${unprocesseddir}/SE_Maps/SE_AP.nii.gz
    #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/SE_Maps/SE_PA.nii.gz -prefix ${unprocesseddir}/SE_Maps/SE_PA.nii.gz
    #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/GRE_FieldMap/GRE_FieldMap.nii.gz -prefix ${unprocesseddir}/GRE_FieldMap/GRE_FieldMap.nii.gz


#fi 
## Reconstruct EPI scans (task and rest)
#if [ -z "$epiRecon" ]; then
    #echo "Skipping EPI reconstruction"
#elif [ $epiRecon = true ]; then
    #echo "Reconstructing Rest EPI Scans..."
    ## make directories for rest scan
    #for ((i=0;i<${#RestEPIs[@]};++i)); do
        ## Count the scans starting from 1 (as opposed to 0)
        #let count=${i}+1

        ## create the unprocessed directory
        #if [ ! -e ${unprocesseddir}/Rest${count} ]; then mkdir ${unprocesseddir}/Rest${count}; fi 
        #echo "Converting Scan ${i} for Rest scan"
        ## reconstruct using mri_convert (freesurfer) for both the epi scan and the sbref
        #mri_convert ${rawdatadir}/${RestEPIs[i]}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/Rest${count}/Rest${count}.nii.gz
        #mri_convert ${rawdatadir}/${RestEPIs_SBRef[i]}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/Rest${count}/Rest${count}_SBRef.nii.gz
        
        ## resample the data to RPI coordinates
        #echo "Resampling Rest images to RPI coordinates"
        #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/Rest${count}/Rest${count}.nii.gz -prefix ${unprocesseddir}/Rest${count}/Rest${count}.nii.gz
        #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/Rest${count}/Rest${count}_SBRef.nii.gz -prefix ${unprocesseddir}/Rest${count}/Rest${count}_SBRef.nii.gz
    #done

    #echo "Reconstructing Task EPI Scans..."
    ## make directories for task scan
    #for ((i=0;i<${#TaskEPIs[@]};++i)); do
        ## Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
        #let count=${i}+1

        ## create the unprocessed directory
        #if [ ! -e ${unprocesseddir}/Task${count} ]; then mkdir ${unprocesseddir}/Task${count}; fi 
        #echo "Converting Scan ${i} for Task scan"
        ## reconstruct using mri_convert (freesurfer) for both the epi scan and the sbref
        #mri_convert ${rawdatadir}/${TaskEPIs[i]}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/Task${count}/Task${count}.nii.gz
        #mri_convert ${rawdatadir}/${TaskEPIs_SBRef[i]}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/Task${count}/Task${count}_SBRef.nii.gz
        
        ## resample the data to RPI coordinates
        #echo "Resampling Task images to RPI coordinates"
        #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/Task${count}/Task${count}.nii.gz -prefix ${unprocesseddir}/Task${count}/Task${count}.nii.gz
        #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/Task${count}/Task${count}_SBRef.nii.gz -prefix ${unprocesseddir}/Task${count}/Task${count}_SBRef.nii.gz
    #done

    #echo "Reconstructing Task Localizer"
    #if [ ! -e ${unprocesseddir}/TaskLoc/ ]; then mkdir ${unprocesseddir}/TaskLoc; fi
    #mri_convert ${rawdatadir}/${TaskLoc}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/TaskLoc/TaskLoc.nii.gz
    #mri_convert ${rawdatadir}/${TaskLoc_SBRef}/*0001.dcm --in_type siemens --out_type nii ${unprocesseddir}/TaskLoc/TaskLoc_SBRef.nii.gz
    #echo "Resampling Task Localizer images to RPI coordinates"
    #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/TaskLoc/TaskLoc.nii.gz -prefix ${unprocesseddir}/TaskLoc/TaskLoc.nii.gz
    #3dresample -overwrite -orient RPI -inset ${unprocesseddir}/TaskLoc/TaskLoc_SBRef.nii.gz -prefix ${unprocesseddir}/TaskLoc/TaskLoc_SBRef.nii.gz



#fi
##########################


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
        --t1="${unprocesseddir}/anat/${subj}_rec-2_T1w.nii.gz" \
        --t2="${unprocesseddir}/anat/${subj}_rec-2_T2w.nii.gz" \
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

#########################
# Start second node of HCP Pipeline: FreeSurferPipeline
if [ -z "$FS" ]; then
    echo "Skipping Freesurfer HCP node"
elif [ $FS = true ]; then
    # limit number of threads of FS
    export OMP_NUM_THREADS=3
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


#########################
# Start fourth node of HCP Pipeline: GenericfMRIVolumeProcessing
#*one block for each of the four tasks: Pracrest, Practask, Testrest, Testtask
if [ -z "$fmriVol" ]; then
    echo "Skipping fMRIVolumeProcessing node"
elif [ $fmriVol = true ]; then
    # need to iterate through each rest scan, and then each task scan
    
    #*note in this version you have to specify --biascorrection method; following parm is taken from the HCP volume example script; previous lab version of the pipeline might have used LEGACY, but SEBASED is better
    BiasCorrection="SEBASED" #NONE, LEGACY, or SEBASED: LEGACY uses the T1w bias field, SEBASED calculates bias field from spin echo images (which requires TOPUP distortion correction)
    

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
	
	## Rest Test
#	echo "Running fMRI Volume processing on Rest Test scan" 
#	fmriname="Rest_Test"
#	fmritcs="${unprocesseddir}/func/${subj}_task-resttest_run-1_bold.nii.gz"
#	fmriscout="${unprocesseddir}/func/${subj}_task-resttest_run-1_sbref.nii.gz"
#	fmap_neg_ap="${unprocesseddir}/fmap/${subj}_dir-AP_run-2_epi.nii.gz"
#	fmap_pos_pa="${unprocesseddir}/fmap/${subj}_dir-PA_run-2_epi.nii.gz"
#	${HCPPipe}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh \
#		--path="${datadir}" \
#		--subject="${subj}" \
#		--fmriname="${fmriname}" \
#		--fmritcs="${fmritcs}" \
#		--fmriscout="${fmriscout}" \
#		--SEPhaseNeg="${fmap_neg_ap}" \
#		--SEPhasePos="${fmap_pos_pa}" \
#		--fmapmag="NONE" \
#		--fmapphase="NONE" \
#		--echospacing="$DwellTime_fMRI" \
#		--echodiff="NONE" \
#		--unwarpdir="${unwarpdir}" \
#		--fmrires="$fmrires" \
#		--dcmethod="TOPUP" \
#		--gdcoeffs="NONE" \
#		--printcom="" \
#		--biascorrection=${BiasCorrection} \
#		--topupconfig="${HCPPIPEDIR_Config}/b02b0.cnf"

    ## Viscomp Task Practice Scans, updated by bpk, 7/30/2018
    #first need to generate number of scans (given aborted runs etc)
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
   
    
fi
#########################

#**Added these fMRI vol nodes to deal with aborted jobs due to pre-emption on e.g. Amarel
#########################
# Start fourth node of HCP Pipeline: GenericfMRIVolumeProcessing
#*one block for each of the four tasks: Pracrest, Practask, Testrest, Testtask
if [ -z "$fmriVol_prac" ]; then
    echo "Skipping fMRIVolumeProcessing_prac node"
elif [ $fmriVol_prac = true ]; then
    # need to iterate through each rest scan, and then each task scan
    
    #*note in this version you have to specify --biascorrection method; following parm is taken from the HCP volume example script; previous lab version of the pipeline might have used LEGACY, but SEBASED is better
    BiasCorrection="SEBASED" #NONE, LEGACY, or SEBASED: LEGACY uses the T1w bias field, SEBASED calculates bias field from spin echo images (which requires TOPUP distortion correction)

    ## Task Practice Scans
    #first need to generate number of scans (given aborted runs etc)
    pushd ${unprocesseddir}/func/
    numPracRuns=`ls *taskpractice*bold.nii.gz | wc -l`
    echo "Prac Task runs for $subj are $numPracRuns"
    popd
	
	#convert string to integer and use that to initiate loop
	fmrivols=($fmriVol_prac_run)
	echo "Start fmriVol from prac run $fmriVol_prac_run"
	
    for ((i=$fmrivols;i<=${numPracRuns};++i)); do
        # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
        echo "Running fMRI Volume processing on Prac Task scan #${i}" 
        
        fmriname="Task_Practice${i}"
        fmritcs="${unprocesseddir}/func/${subj}_task-taskpractice_run-${i}_bold.nii.gz"
		fmriscout="${unprocesseddir}/func/${subj}_task-taskpractice_run-${i}_sbref.nii.gz"
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
   
fi
#########################

#########################
if [ -z "$fmriVol_test" ]; then
    echo "Skipping fMRIVolumeProcessing_test node"
elif [ $fmriVol_test = true ]; then    
    #*note in this version you have to specify --biascorrection method; following parm is taken from the HCP volume example script; previous lab version of the pipeline might have used LEGACY, but SEBASED is better
    BiasCorrection="SEBASED" #NONE, LEGACY, or SEBASED: LEGACY uses the T1w bias field, SEBASED calculates bias field from spin echo images (which requires TOPUP distortion correction)
    
    ## Task Test Scans
    #first need to generate number of scans (given aborted runs etc)
    pushd ${unprocesseddir}/func/
    numTestRuns=`ls *tasktest*bold.nii.gz | wc -l`
    echo "Test Task runs for $subj are $numTestRuns"
    popd
    
    #convert string to integer and use that to initiate loop
	fmrivols=($fmriVol_test_run)
	echo "Start fmriVol from test run $fmriVol_test_run"
	
    for ((i=$fmrivols;i<=${numTestRuns};++i)); do
        # Count the scans starting from 1 (as opposed to 0). This is for the numbering for the nifti filename.
        echo "Running fMRI Volume processing on Test Task scan #${i}" 
        
        fmriname="Task_Test${i}"
        fmritcs="${unprocesseddir}/func/${subj}_task-tasktest_run-${i}_bold.nii.gz"
		fmriscout="${unprocesseddir}/func/${subj}_task-tasktest_run-${i}_sbref.nii.gz"
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

#########################
# Start fifth node of HCP Pipeline: GenericfMRISurfaceProcessing
#*one block for each of the four tasks: Pracrest, Practask, Testrest, Testtask
if [ -z "$fmriSurf" ]; then
    echo "Skipping fMRI Surface Processing node"
elif [ $fmriSurf = true ]; then

    ## Rest, by bpk, 8/6/2018
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
	
	## Rest Test
    # set input name
   # fmriname="Rest_Test"
   # ${HCPPipe}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh \
   #		--path="${datadir}" \
   #     --subject="${subj}" \
   #     --fmriname="${fmriname}" \
   #     --fmrires="$fmrires" \
   #     --lowresmesh="$LowResMeshes" \
   #     --smoothingFWHM=$SmoothingFWHM \
   #    --grayordinatesres="$GrayordinatesResolutions" \
   #     --regname=$RegName

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
        
    
fi
#########################


#***Add new modules
#########################


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
    
    ## Rest Test, commented out by bpk, 8/2/18
    # set input name
    #fmriname="Rest_Test"
    #fmri_input=${subjdir}/${subj}/MNINonLinear/Results/${fmriname}/${fmriname}.nii.gz
    #${FSL_FIXDIR}/hcp_fix ${fmri_input} ${bandpass}
    
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
	
	#generate run names for Task_Retino and Task_Viscomp (these vary per subject due to aborted scans)
    pushd ${unprocesseddir}/func/
    numViscompRuns=`ls *viscomp*bold.nii.gz | wc -l`
    echo " Task runs for $subj are $numViscompRuns"
    popd
    
    pushd ${unprocesseddir}/func/
    numRetinoRuns=`ls *retino*bold.nii.gz | wc -l`
    echo "Test Task runs for $subj are $numRetinoRuns"
    popd
    
	run_names="Task_Viscomp1" #init string
    for ((i=2;i<=${numViscompRuns};++i)); do
		run_names="${run_names}@Task_Viscomp${i}"
    done
	
	run_names="${run_names}@Task_Retino1"
	for ((i=2;i<=${numRetinoRuns};++i)); do
		run_names="${run_names}@Task_Retino${i}"
    done
	
	#*@ delimited fMRIName strings corresponding to maps that will have ICA+FIX reapplied to them (could be either rfMRI or tfMRI). If none are to be used, specify "NONE".
	#rfMRINames="rfMRI_REST1_LR rfMRI_REST1_RL rfMRI_REST2_LR rfMRI_REST2_RL" #Space delimited list or NONE
	#rfMRINames="Rest_Practice@Rest_Test@Task_Practice@Task_Test"
	rfMRINames="NONE"
	#*@ delimited fMRIName strings corresponding to maps that will not have ICA+FIX reapplied to them (likely not to be used in the future as ICA+FIX will be recommended for all fMRI data) If none are to be used, specify "NONE".
	#tfMRINames="tfMRI_WM_LR tfMRI_WM_RL tfMRI_GAMBLING_LR tfMRI_GAMBLING_RL tfMRI_MOTOR_LR tfMRI_MOTOR_RL tfMRI_LANGUAGE_LR tfMRI_LANGUAGE_RL tfMRI_SOCIAL_LR tfMRI_SOCIAL_RL tfMRI_RELATIONAL_LR tfMRI_RELATIONAL_RL tfMRI_EMOTION_LR tfMRI_EMOTION_RL" #Space delimited list or NONE
	#*not sure if run numbers need to be added for task scans?
	tfMRINames="Task_Rest@${run_names}"
	echo "tfMRINames for this subject are $tfMRINames"
	
	#Should equal previous grayordiantes smoothing in fMRISurface (because we are resampling from unsmoothed native mesh timeseries - already set
	#SmoothingFWHM="2" 
	#HighPass=${bandpass}
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
   


#**Do not run - everything from here on out is handled by Mike's matlab GLM scripts (i.e. func runs concatenation, mask extraction etc)
#########################
# Create masks
if [ -z "$createMasks" ]; then
    echo "Skipping mask creation"
elif [ "$createMasks" == "true" ]; then
    
    echo "Creating gray, white, ventricle, whole brain masks for subject ${subj}..."

    # HCP standard to parcel out white v gray v ventricle matter
    segparc=${subjdir}/MNINonLinear/wmparc.nii.gz

    # Change to subjmaskdir
    pushd $subjmaskdir

    
    
    ###############################
    ### Create whole brain masks
    echo "Creating whole brain mask for subject ${subj}..."
    3dcalc -overwrite -a $segparc -expr 'ispositive(a)' -prefix ${subj}_wholebrainmask.nii.gz
    # Resample to functional space (func in MNI space)
    3dresample -overwrite -master ${subjdir}/MNINonLinear/Results/Rest_Practice/Rest_Practice.nii.gz -inset ${subj}_wholebrainmask.nii.gz -prefix ${subj}_wholebrainmask_func.nii.gz
    # Dilate mask by 1 functional voxel (just in case the resampled anatomical mask is off by a bit)
    3dLocalstat -overwrite -nbhd 'SPHERE(-1)' -stat 'max' -prefix ${subj}_wholebrainmask_func_dil1vox.nii.gz ${subj}_wholebrainmask_func.nii.gz
    
   

    ###############################
    ### Create gray matter masks
    echo "Creating gray matter masks for subject ${subj}..." 
    # Indicate the mask value set for wmparc.nii.gz
    # Gray matter mask set
    maskValSet="8 9 10 11 12 13 16 17 18 19 20 26 27 28 47 48 49 50 51 52 53 54 55 56 58 59 60 96 97 1000 1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1017 1018 1019 1020 1021 1022 1023 1024 1025 1026 1027 1028 1029 1030 1031 1032 1033 1034 1035 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035"

    # Add segments to mask
    maskNum=1
    for maskval in $maskValSet
    do
	if [ ${maskNum} = 1 ]; then
            3dcalc -a $segparc -expr "equals(a,${maskval})" -prefix ${subj}mask_temp.nii.gz -overwrite
        else
            3dcalc -a $segparc -b ${subj}mask_temp.nii.gz -expr "equals(a,${maskval})+b" -prefix ${subj}mask_temp.nii.gz -overwrite
        fi
	let maskNum++
    done
    #Make mask binary
    3dcalc -a ${subj}mask_temp.nii.gz -expr 'ispositive(a)' -prefix ${subj}_gmMask.nii.gz -overwrite
    #Resample to functional space
    3dresample -overwrite -master ${subjdir}/MNINonLinear/Results/Rest_Practice/Rest_Practice.nii.gz -inset ${subj}_gmMask.nii.gz -prefix ${subj}_gmMask_func.nii.gz
    #Dilate mask by 1 functional voxel (just in case the resampled anatomical mask is off by a bit)
    3dLocalstat -overwrite -nbhd 'SPHERE(-1)' -stat 'max' -prefix ${subj}_gmMask_func_dil1vox.nii.gz ${subj}_gmMask_func.nii.gz
    
    rm -f ${subj}mask_temp.nii.gz
       
   
    
    ###############################
    ### Create white matter masks
    echo "Creating white matter masks for subject ${subj}..."

    # Indicate the mask value set for wmparc.nii.gz
    # White matter mask set
    maskValSet="250 251 252 253 254 255 3000 3001 3002 3003 3004 3005 3006 3007 3008 3009 3010 3011 3012 3013 3014 3015 3016 3017 3018 3019 3020 3021 3022 3023 3024 3025 3026 3027 3028 3029 3030 3031 3032 3033 3034 3035 4000 4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033 4034 4035 5001 5002"

    # Add segments to mask
    maskNum=1
    for maskval in $maskValSet
    do
	if [ ${maskNum} = 1 ]; then
            3dcalc -a $segparc -expr "equals(a,${maskval})" -prefix ${subj}mask_temp.nii.gz -overwrite
        else
            3dcalc -a $segparc -b ${subj}mask_temp.nii.gz -expr "equals(a,${maskval})+b" -prefix ${subj}mask_temp.nii.gz -overwrite
        fi
	let maskNum++
    done
    #Make mask binary
    3dcalc -a ${subj}mask_temp.nii.gz -expr 'ispositive(a)' -prefix ${subj}_wmMask.nii.gz -overwrite
    #Resample to functional space
    3dresample -overwrite -master ${subjdir}/MNINonLinear/Results/Rest_Practice/Rest_Practice.nii.gz -inset ${subj}_wmMask.nii.gz -prefix ${subj}_wmMask_func.nii.gz
    #Subtract graymatter mask from white matter mask (avoiding negative #s)
    3dcalc -a ${subj}_wmMask_func.nii.gz -b ${subj}_gmMask_func_dil1vox.nii.gz -expr 'step(a-b)' -prefix ${subj}_wmMask_func_eroded.nii.gz -overwrite
    
    rm -f ${subj}mask_temp.nii.gz
          

    
    ###############################
    ### Create ventricle masks
    echo "Creating ventricle matter masks for subject ${subj}..."

    # Indicate the mask value set for wmparc.nii.gz
    # Ventricle mask set
    maskValSet="4 43 14 15"

    # Add segments to mask
    maskNum=1
    for maskval in $maskValSet
    do
	if [ ${maskNum} = 1 ]; then
            3dcalc -a $segparc -expr "equals(a,${maskval})" -prefix ${subj}mask_temp.nii.gz -overwrite
        else
            3dcalc -a $segparc -b ${subj}mask_temp.nii.gz -expr "equals(a,${maskval})+b" -prefix ${subj}mask_temp.nii.gz -overwrite
        fi
	let maskNum++
    done
    #Make mask binary
    3dcalc -a ${subj}mask_temp.nii.gz -expr 'ispositive(a)' -prefix ${subj}_ventricles.nii.gz -overwrite
    #Resample to functional space
    3dresample -overwrite -master ${subjdir}/MNINonLinear/Results/Rest_Practice/Rest_Practice.nii.gz -inset ${subj}_ventricles.nii.gz -prefix ${subj}_ventricles_func.nii.gz
    #Subtract graymatter mask from ventricles (avoiding negative #s)
    3dcalc -a ${subj}_ventricles_func.nii.gz -b ${subj}_gmMask_func_dil1vox.nii.gz -expr 'step(a-b)' -prefix ${subj}_ventricles_func_eroded.nii.gz -overwrite
    rm -f ${subjNum}mask_temp.nii.gz
    
    rm -f ${subj}mask_temp.nii.gz
          
    popd

fi 



#########################
# Concatenate Task Runs - *not sure if this is completely necessary (will have to check Mike's GLM script to make sure), but run it just for consistency; **modify to write correct run numbers for each subject
#*will need to factor in aborted runs
if [ -z "$concatenateRuns" ]; then
    echo "Skipping run concatenation"
elif [ "$concatenateRuns" == "true" ]; then
    
    pushd ${analysisdir}

    echo "Concatenating the 8 task runs for subject ${subj}..."
    runs="1 2 3 4 5 6 7 8" 
    runList=""
    for run in $runs
    do
        runList="${runList} ${subjdir}/MNINonLinear/Results/Task${run}/Task${run}.nii.gz"
    done

    3dTcat -prefix Task_allruns.nii.gz ${runList}

    popd
fi
  
#########################
# Extract Timeseries - **modify so that timeseries are written for each of the 4 scans; also takes concat timeseries as input
if [ -z "$tsExtract" ]; then
    echo "Skipping Timeseries Extraction"
elif [ "$tsExtract" == "true" ]; then

    pushd ${analysisdir}

    echo "Extracting timeseries from Rest1.nii.gz and Task_allruns.nii.gz using white matter masks, ventricles, and whole brain masks for subject ${subj}..."
    rest=${subjdir}/MNINonLinear/Results/Rest1/Rest1.nii.gz
    task=${analysisdir}/Task_allruns.nii.gz
    
    # Extract Rest TS using wm and ventricle masks
    3dmaskave -quiet -mask ${subjmaskdir}/${subj}_wmMask_func_eroded.nii.gz ${rest} > ${subj}_WM_timeseries_rest.1D
    3dmaskave -quiet -mask ${subjmaskdir}/${subj}_ventricles_func_eroded.nii.gz ${rest} > ${subj}_ventricles_timeseries_rest.1D
    3dmaskave -quiet -mask ${subjmaskdir}/${subj}_wholebrainmask_func_dil1vox.nii.gz ${rest} > ${subj}_wholebrainsignal_timeseries_rest.1D

    # Extract Task TS using wm and ventricle masks
    3dmaskave -quiet -mask ${subjmaskdir}/${subj}_wmMask_func_eroded.nii.gz ${task} > ${subj}_WM_timeseries_task.1D
    3dmaskave -quiet -mask ${subjmaskdir}/${subj}_ventricles_func_eroded.nii.gz ${task} > ${subj}_ventricles_timeseries_task.1D
    3dmaskave -quiet -mask ${subjmaskdir}/${subj}_wholebrainmask_func_dil1vox.nii.gz ${task} > ${subj}_wholebrainsignal_timeseries_task.1D
    
    popd
fi
    
    
#*Exclude this
######################### DOES THIS NEED TO BE DONE? OR DOES HCP ALREADY COVER THIS?
## Resting-state Motion Correction
#if [ -z "$restMotionCorrectoin" ]; then
#    echo "Skipping Timeseries Extraction"
#elif [ "$tsExtract" == "true" ]; then
#
#fi
#
#
#########################
# Task-based analysis!
# Start Run concatenation and GLM analysis
execute=0
if [ $execute = true ]; then
    # First need to concat all 8 Task Runs together
    # Create run list
    runList=""
    concatString="1D:"
    TRCount=0
    for ((i=1;i<=${#TaskEPIs[@]};++i)); do
        runList="${runList} ${subjdir}/MNINonLinear/Results/Task${i}/Task${i}.nii.gz"
        concatString="${concatString} ${TRCount}"
        TRCount=$(expr $TRCount + $numTRsPerTaskRun)
    done

    echo "-Concatenating task runs-"
    echo "Run list: ${runList}"
    echo "Concatenation string (onset times of each run): ${concatString}"
    # create Analysis directory
    if [ ! -e ${subjdir}/Analysis ]; then mkdir ${subjdir}/Analysis; fi 
    #3dTcat -prefix ${subjdir}/Analysis/Task_allruns ${runList}

    echo "-Running GLM-"
   
    # First resample freesurfer brain mask to functional space (From HCP Pipelines)
    3dresample -master ${subjdir}/Analysis/Task_allruns+tlrc -prefix ${subjdir}/Analysis/brainmask_fs_func.nii.gz -inset ${subjdir}/MNINonLinear/brainmask_fs.nii.gz -overwrite

    3dDeconvolve \
        -input ${subjdir}/Analysis/Task_allruns+tlrc \
        -concat "$concatString" \
        -mask ${subjdir}/Analysis/brainmask_fs_func.nii.gz \
        -polort A \
        -num_stimts 1 \
        -stim_times 1 ${subjdir}/timingfiles/stime_013_stimfile_IndivRITL_PilotGLM_EV1_task.1D.01.1D 'BLOCK(1,1)' -stim_label 1 Task \
        -errts ${subjdir}/Analysis/residual_error_series \
        -fout -tout \
        -xsave -x1D ${subjdir}/Analysis/xmat_rall.x1D -xjpeg ${subjdir}/Analysis/xmat_rall.jpg \
        -jobs 8 -float -overwrite \
        -bucket ${subjdir}/Analysis/pilotGLM_outbucket -cbucket ${subjdir}/Analysis/pilotGLM_cbucket
fi

