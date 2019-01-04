 function [output] = Step2b_GLM_PostHCPpreproc_bpk(SUBJECTLIST)

%% Models
% Model1: Block model--All frag (including encoding) versus all illus (including encoding)
% Model2: Encode trials (instructions) or stimulus trials (2 regressors)
% Model3: Encode (frag or illus), and stimtulus (frag/illus) (4 regressors)
% Model4: Same as Model3, except that RT decides "duration" of stimulus
% Model5: New: Same as Model3, except that RT is used as additional regressor (6 regressors)
% Model6: New: Same as model 3, but does not use regressors for motion

% This function runs preprocessing and GLM analysis, after the HCP minimal preproccessing pipeline has been run.
% It is designed to parcellate a dataset into a set of regions, which are then preprocessed.
% Ravis Script
% Preprocessing steps included:
% 1) Parcellate a dense CIFTI file into a set of N time series, where N is the number of parcels/regions
% 2) Prepare nuisance regressors for removal. This includes spatial mask definition (e.g., white matter)
%    and extraction from volume fMRI data, as well as processing of motion parameters.
% 3) Preparation of task regressors, when task runs are present. This includes (with custom scripts)
%    conversion of behavioral timing data to a common format, convolving with a hemodynamic response function, 
%    and conversion to a GLM design matrix.
% 4) Rest fMRI nuisance regression for functional connectivity analyses, if rest data are present.
% 5) Task fMRI GLM along with nuisance regression (if task data are present), for functional connectivity and/or task activation analyses.
% 6) Temporal filtering of time series (optional).
%
% Note: Frequently customized variables are IN CAPS throughout the script

%Script author:
%Michael W. Cole, mwcole@mwcole.net, http://www.colelab.org

%Script version: 1.5
%Date: October 25, 2016

% Modified by Ravi Mil
%
%Version history:

%RM Edits for HCP pipeline v2 (with MSMsulc/MSMall)
%1. Bandpass filter (2000s) was not applied to tfmri files in last script of
%pipeline, so current data has not been detrended - run regressors should achieve this?
%Regardless seems like motion regs are not affected by filter (nb: _dt.txt in the preproc output are the detrended motion
%regs; use the standard .txt regressors as mentioned in Mike's script
%below). Mike's filtering only applies to rest data - apply this as per
%parms below? *DECIDE ON: filtering and motion scrubbing (for task and
%rest)? Leave both off for task data (interpolation is used for scrubbed TRs); apply both to rest data?

%***=filtering OFF for task+rest (assume GLM linear detrend will work);motion regressor filtering ON for task+rest
%=motion scrubbing OFF for task, ON for rest

%OVERALL PLAN 
% 1. Use ciftiopen to open a dscalar.nii as a template structure: cii = ciftiopen('path/to/file','path/to/wb_command'); CIFTIdata = cii.cdata;
% 2. Run Mike's GLM script (see NOTES immediately above for parm choices):
% 2a. Parcellate the dtseries file output from HCP preproc using the Glasser parcels (stored in scripts dir). This is done separately for L/R surface hemispheres, and the result concatenated (1-180=L; 181-360=R).
% 2b. Run the GLM - *make sure that timingFiles are read in correctly
% 3. Group analysis
    % Compute contrast images for each subject if necessary (for Model2 only)
    % Average betas (or contrast images) for each subject.
    % Convert the group tstat to zstat (using andybrainblog code): norminv(tcdf(i1,' num2str(dof) '),0,1), where i1=tstat, dof=numsubs-1
    % In AnalysisOutput, store activation maps for zstats thresholded/binarized by significance across following thresholds: p < .05, FDR p < .05, raw.
    % *Also store the validation result metric:
        % Plot num voxels with significant activation (based on various p thresholds) = 12 GLMs (model1a,1b,1c,2a,2b,2c x MSMall/MSMsulc) x separately for the 2 prac/test sessions.
        % Plot peak zstats = 12 GLMs x separately for the 2 prac/test sessions.
    % Output format:
    % AnalysisOutput.(model).(MSMversion).(session).(threshold).beta/zstat/num_sig_voxels/peak_zstat
% 4. Write new cifti files using AnalysisOutput: newcii = cii; newcii.cdata = AnalysisOutput.(session).(threshold); ciftisave(newcii,'path/to/newfile','path/to/wb_command').
    % Output should be a vector of zstat amplitudes - one for each 32k l/r hemi grayordinate (rather than one for each 360 glasser region). This means that I will need the dlabel file coding for which 'grayordinate' corresponds to which Glasser region, to enter a single region amplitude value for all grayordinates for the 360 regions.
    % Might have to use ciftisavereset as the data matrices will have a different number of maps/columns from what you started with: ciftisavereset(newcii,'path/to/newfile','path/to/wb_command');
% 5. Use workbench to visualize - use the Conte 32k atlas as underlay, and the zstat images as the overlay.



% trials within block

%timeStart=tic;
%% Parameters to customize for your analysis
addpath('/projects/AnalysisTools/')
addpath('/projects/AnalysisTools/gifti-1.6/')
addpath('/projects/AnalysisTools/ReffuncConverter/')
%Dataset-specific paths
addpath('/projects3/NeuralMech/docs/docs/scripts/')

%%Basic processing parameters
if nargin<1
     SUBJECTLIST={'sub-C05','sub-C13','sub-C14','sub-C15','sub-C22','sub-C25','sub-C26','sub-C27', 'sub-C29','sub-C32', 'sub-C33'};
   % SUBJECTLIST={'sub-C05'};
   %SUBJECTLIST={'sub-C05','sub-C13','sub-C14','sub-C15','sub-C22'};
end
numSubjs=length(SUBJECTLIST);
TR_INSECONDS=0.785;
FRAMESTOSKIP=0; % should ideally skip first 5 TRs for Rest only

%Basic data parameters - 
%RUNNAMES is determined by names given to func runs by HCP preproc...

func_info=struct;
% Be sure to put "rest" run first below:
%testinfo.RUNNAMES = {'Task_Rest','Task_Retino1','Task_Retino2','Task_Retino3','Task_Viscomp1','Task_TViscomp2','Task_Viscomp3','Task_Viscomp4'};
func_info.RUNNAMES = {'Task_Rest','Task_Viscomp1','Task_Viscomp2','Task_Viscomp3','Task_Viscomp4'};
func_info.numRuns=length(func_info.RUNNAMES);
func_info.RESTRUNS=1;
func_info.TASKRUNS=2:5;
%func_info.RUNLENGTHS = [765, 306, 306, 306, 281, 281, 281, 281];
func_info.RUNLENGTHS = [765, 281, 281, 281, 281];

%parcel info
L_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
R_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
parcellationName='Glasser2016';
NUMPARCELS=360;

%Data processing flags
GSR=0;      %GSR = 0 if no GSR, 1 if you want to include GSR
NPROC=4;    %Number of processors to use when computing GLMs
FDTHRESH=0.30;  %Framewise displacement; The threshold (in millimeters) for flagging in-scanner movement; 0.3 was used in Doug's 2018 paper
TEMPORALFILTER=0;   %Flag indicating if data should be temporally filtered
HIGHPASS_HZ=0.008;
LOWPASS_HZ=0.09;
IMPLEMENT_MOTIONSCRUBBING=2;     %Set to 1 for yes, 0 for no; **2=should 1 for rest, 0 for task; MODIFY CODE BELOW!

%Directories
BASEDIR='/projects3/NeuralMech/';
timingfileDir=[BASEDIR 'data/timingFiles/']; %these files give regressors
datadir=[BASEDIR 'data/preprocessed/'];
outputdatadir=[BASEDIR '/data/results/GLM/'];
if ~exist(outputdatadir, 'dir'); mkdir(outputdatadir); end

%**Set up Model names (should help input StimFiles + name output struct
% stim_model_input={'Model1a_1Taskreg_varblock','Model1b_1Taskreg_consepochlength','Model1c_1Taskreg_varepochRT',...
%    'Model2a_64Taskreg_varblock','Model2b_64Taskreg_consepochlength','Model2c_64Taskreg_varepochRT'};
stim_model_input={'model1','model2','model3', 'model4', 'model5', 'model6'};

%set input suffix which will be chosen based on loop below (between MSMsulc
%and MSMsulc_MSMall
%MSMsulc output=$sub/MNINonLinear/Results/$RUN/$RUN_Atlas.dtseries.nii
%MSMall output=$sub/MNINonLinear/Results/$RUN/$RUN_Atlas_MSMAll_InitalReg_2_d40_WRN.dtseries.nii
MSM_input={'_Atlas.dtseries.nii','_Atlas_MSMAll_InitalReg_2_d40_WRN.dtseries.nii'};
MSM_struct_name={'sulc','msmAll'};

%**determines output filename - set based on 12 GLM variations (with
%additional suffix set by Prac/Test)
ANALYSISNAME='GLMs_Viscomp';

%% Iterate through subjects 
%Start loops: model->MSM version
GLMOutput=struct;
execute=1; % set to 0 if GLM has already been run

if execute==1
    for modeli=1:3%length(stim_model_input)
        for MSM_version=1:length(MSM_input)
            %set names (eg for output name)
            model_name=stim_model_input{modeli};
            MSM_name=MSM_struct_name{MSM_version};

            %start subject loop
            output=[];
            output.SUBJECTLIST=SUBJECTLIST;
            output.RUNNAMES=func_info.RUNNAMES;
            output.RUNLENGTHS=func_info.RUNLENGTHS;
            output.RESTRUNS=func_info.RESTRUNS;
            output.TASKRUNS=func_info.TASKRUNS;
            output.parcellationName=parcellationName;

            for subjIndex=1:numSubjs
                subjName = SUBJECTLIST{subjIndex};
                disp(['Processing subject ' subjName]);

                %% Modify parms for subjects with missing/aborted runs                            
%                   if strcmp(subjNum,'sub-6')
%                         func_info.RUNNAMES = {'Rest_Test','Task_Test1','Task_Test2','Task_Test3','Task_Test4','Task_Test5','Task_Test6','Task_Test7','Task_Test8'};
%                         func_info.numRuns=length(func_info.RUNNAMES);
%                         func_info.RESTRUNS=1;
%                         func_info.TASKRUNS=2:9;
%                         func_info.RUNLENGTHS = [1147, 521, 521, 521, 521, 521, 521, 521, 521];
%                   end
                
                tseriesMatSubj=zeros(NUMPARCELS,max(func_info.RUNLENGTHS),func_info.numRuns); % 3 dimensional matrix
                
                %Directories
                subjDir=[datadir '/' subjName '/'];
                if ~exist(subjDir, 'dir'); mkdir(subjDir); end
                subjTemporaryAnalysisDir=[subjDir ANALYSISNAME '/' model_name '_' MSM_name '/'];
                if ~exist(subjTemporaryAnalysisDir, 'dir'); mkdir(subjTemporaryAnalysisDir); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% DOWNSAMPLING GRAYORDINATE DATA TO PARCELS (GLASSER REGIONS) %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Downsampling grayordinate data to parcels') % creates 3D matrix tSeriesMatSubj, which is plugged into rest and task GLMs below

                %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
                RERUN_EXTRACTPARCELS=1;

                runCount=1;
                for runName=func_info.RUNNAMES
                    thisRunName=runName{1};
                    % These are the parcellated filenames
                    L_parcelTSFilename=[subjTemporaryAnalysisDir '/' thisRunName '_Atlas.L.' parcellationName 'Parcels.32k_fs_LR.ptseries.nii'];
                    R_parcelTSFilename=[subjTemporaryAnalysisDir '/' thisRunName '_Atlas.R.' parcellationName 'Parcels.32k_fs_LR.ptseries.nii'];

                    %**set subj_MSM based on MSM_version
                    subj_MSM=MSM_input{MSM_version};

                    if RERUN_EXTRACTPARCELS == 1
                        subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
                        inputFile=[subjRunDir thisRunName subj_MSM]; % Below, the cifti files are from HCP
                        eval(['!wb_command -cifti-parcellate ' inputFile ' ' L_parcelCIFTIFile ' COLUMN ' L_parcelTSFilename ' -method MEAN'])
                        eval(['!wb_command -cifti-parcellate ' inputFile ' ' R_parcelCIFTIFile ' COLUMN ' R_parcelTSFilename ' -method MEAN'])
                    end

                    %Load newly generated parcellated data for this run
                    L_dat = ciftiopen(L_parcelTSFilename,'wb_command');
                    R_dat = ciftiopen(R_parcelTSFilename,'wb_command');
                    if size(L_dat.cdata,2)>func_info.RUNLENGTHS(runCount)
                        disp(['WARNING: More TRs for this run than expected. Subject: ' subjName ', Run: ' num2str(runCount)])
                    elseif size(L_dat.cdata,2)<func_info.RUNLENGTHS(runCount)
                        disp(['WARNING: Fewer TRs for this run than expected. Subject: ' subjName ', Run: ' num2str(runCount)])
                    end
                    % Each row = parcel, each column = volume, each slice = run
                    tseriesMatSubj(1:(NUMPARCELS/2),1:func_info.RUNLENGTHS(runCount),runCount)=L_dat.cdata(:,1:func_info.RUNLENGTHS(runCount));
                    tseriesMatSubj((NUMPARCELS/2+1):end,1:func_info.RUNLENGTHS(runCount),runCount)=R_dat.cdata(:,1:func_info.RUNLENGTHS(runCount));
                    runCount=runCount+1;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% PREPARE NUISANCE REGRESSORS %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Preparing nuisance regressors')
                %Using Freesurfer aparc+aseg masks

                %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
                RERUN_PREPNUISANCEREG=1;
                savedNuisRegfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_nuisanceTSVars.mat'];

                if or(RERUN_PREPNUISANCEREG == 1, ~exist(savedNuisRegfile, 'file'))
                    subjMaskDir=[subjDir '/masks/'];
                    if ~exist(subjMaskDir, 'dir'); mkdir(subjMaskDir); end

                    %Resample Freesurfer segmented mask into functional space (using nearest neighbor interpolation); uses AFNI
                    subjMNINonLinearDir=[subjDir '/MNINonLinear/'];
                    exampleFunctionalVolFile=[subjMNINonLinearDir '/Results/' func_info.RUNNAMES{1} '/' func_info.RUNNAMES{1} '.nii.gz'];
                    eval(['!3dresample -overwrite -rmode NN -master ' exampleFunctionalVolFile ' -inset ' subjMNINonLinearDir 'aparc+aseg.nii.gz -prefix ' subjMNINonLinearDir 'aparc+aseg_resampFunc.nii.gz']);

                    %Load Freesurfer segmented mask
                    aparc_aseg=load_nifti([subjMNINonLinearDir '/aparc+aseg_resampFunc.nii.gz']);

                    %Create gray matter mask
                    maskValSet_graymatter=[8 9 10 11 12 13 16 17 18 19 20 26 27 28 47 48 49 50 51 52 53 54 55 56 58 59 60 96 97 1000:1035 2000:2035];
                    grayMatterMask=ismember(aparc_aseg.vol,maskValSet_graymatter);

                    %Create white matter mask
                    maskValSet_whitematter=[2 7 41 46];
                    whiteMatterMask=ismember(aparc_aseg.vol,maskValSet_whitematter);
                    %Erode white matter mask by 2 voxels
                    whiteMatterMask_eroded=imerode(whiteMatterMask,strel(ones(2,2,2)));

                    %Create ventricle mask
                    maskValSet_ventricles=[4 43 14 15];
                    ventricleMask=ismember(aparc_aseg.vol,maskValSet_ventricles);
                    %Erode ventricle matter mask by 2 voxels
                    ventricleMask_eroded=imerode(ventricleMask,strel(ones(2,2,2)));

                    %Create whole brain mask
                    wholebrainMask=aparc_aseg.vol>0;

                    %Load in nuisance time series for each run
                    nuisanceTS_whitematter=zeros(max(func_info.RUNLENGTHS),func_info.numRuns);
                    nuisanceTS_ventricles=zeros(max(func_info.RUNLENGTHS),func_info.numRuns);
                    nuisanceTS_wholebrain=zeros(max(func_info.RUNLENGTHS),func_info.numRuns);
                    nuisanceTS_motion=zeros(12,max(func_info.RUNLENGTHS),func_info.numRuns);
                    FD_motion=zeros(max(func_info.RUNLENGTHS),func_info.numRuns);
                    temporalMask=ones(max(func_info.RUNLENGTHS),func_info.numRuns);
                    numFramesCensored=zeros(func_info.numRuns,1);

                    runCount=1;
                    for runName=func_info.RUNNAMES
                        thisRunName=runName{1};
                        subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
                        %inputFile=[subjRunDir thisRunName subj_MSM];
                        inputFile=[subjRunDir thisRunName '.nii.gz'];

                        %Load data
                        runData=load_nifti(inputFile);
                        runData2D=reshape(runData.vol,size(runData.vol,1)*size(runData.vol,2)*size(runData.vol,3),size(runData.vol,4));
                        if size(runData2D,2)>func_info.RUNLENGTHS(runCount)
                            disp(['WARNING: More TRs for this run than expected. Subject: ' subjName ', Run: ' thisRunName ', Run number: ' num2str(runCount)])
                            runData2D=runData2D(:,1:func_info.RUNLENGTHS(runCount));
                        elseif size(runData2D,2)<func_info.RUNLENGTHS(runCount)
                            disp(['WARNING: Fewer TRs for this run than expected. Subject: ' subjName ', Run: ' thisRunName ', Run number: ' num2str(runCount)])
                            runData2D=runData2D(:,1:func_info.RUNLENGTHS(runCount));
                        end

                        whiteMatterMask_eroded_1D=reshape(whiteMatterMask_eroded,size(whiteMatterMask_eroded,1)*size(whiteMatterMask_eroded,2)*size(whiteMatterMask_eroded,3),1);
                        nuisanceTS_whitematter(1:func_info.RUNLENGTHS(runCount),runCount)=mean(runData2D(whiteMatterMask_eroded_1D,:),1);

                        ventricleMask_eroded_1D=reshape(ventricleMask_eroded,size(ventricleMask_eroded,1)*size(ventricleMask_eroded,2)*size(ventricleMask_eroded,3),1);
                        nuisanceTS_ventricles(1:func_info.RUNLENGTHS(runCount),runCount)=mean(runData2D(ventricleMask_eroded_1D,:),1);

                        wholebrainMask_1D=reshape(wholebrainMask,size(wholebrainMask,1)*size(wholebrainMask,2)*size(wholebrainMask,3),1);
                        nuisanceTS_wholebrain(1:func_info.RUNLENGTHS(runCount),runCount)=mean(runData2D(wholebrainMask_1D,:),1);

                        %Note: derivatives are already included in motion time series
                        motionvals=importdata([subjRunDir 'Movement_Regressors.txt'])';
                        nuisanceTS_motion(:,1:func_info.RUNLENGTHS(runCount),runCount)=motionvals(:,1:func_info.RUNLENGTHS(runCount));

                        %Skip first FRAMESTOSKIP frames
                        temporalMask(1:FRAMESTOSKIP,runCount)=0;

                        %Calculate framewise displacement (FD) according to Power et al. (2012)
                        %Briefly: The sum of the absolute values of the translational and rotational displacements over all frames (in mm)
                        %Note: HCP's minimal preprocessing pipeline uses the following ordering of the motion parameters (see https://github.com/Washington-University/Pipelines/blob/master/global/scripts/mcflirt_acc.sh):
                        %trans x, trans y, trans z, rot x, rot y, rot z [rotations in degrees], then derivatives of those 6 (for 12 total)
                        motionTS_dt=[zeros(size(nuisanceTS_motion,1),1) diff(squeeze(nuisanceTS_motion(:,1:func_info.RUNLENGTHS(runCount),runCount))')'];
                        assumedRadius=50;
                        rot_x=(2*assumedRadius*pi/360)*motionTS_dt(4,:);
                        rot_y=(2*assumedRadius*pi/360)*motionTS_dt(5,:);
                        rot_z=(2*assumedRadius*pi/360)*motionTS_dt(6,:);
                        FD_motion(1:func_info.RUNLENGTHS(runCount),runCount)=abs(motionTS_dt(1,:))+abs(motionTS_dt(2,:))+abs(motionTS_dt(3,:))+abs(rot_x)+abs(rot_y)+abs(rot_z);
                        %Apply temporal filtering to FD, to reduce the effect of respiration on FD measure. Based on Siegel et al. (2016) [Siegel JS, Mitra A, Laumann TO, Seitzman BA, Raichle M, Corbetta M, Snyder AZ (2016) ?Data Quality Influences Observed Links Between Functional Connectivity and Behavior?. Cereb Cortex. 1?11.http://doi.org/10.1093/cercor/bhw253]
                        %Create temporal filter
                        lopasscutoff=0.3/(0.5/TR_INSECONDS); % lowpass filter of 0.3 Hz
                        %Using filter order of 1
                        filtorder=1;
                        [butta, buttb]=butter(filtorder,lopasscutoff);
                        %Apply temporal filter
                        filteredFD=filtfilt(butta,buttb,FD_motion(1:func_info.RUNLENGTHS(runCount),runCount));
                        FD_motion(1:func_info.RUNLENGTHS(runCount),runCount)=filteredFD;

                        %Implement motion scrubbing/censoring in the temporal mask
                        %**set IMPLEMENT_MOTIONSCRUBBING for REST runs only!
                        if IMPLEMENT_MOTIONSCRUBBING==2
                            if runCount==1 %rest is always first
                                threshedMotion=FD_motion(1:func_info.RUNLENGTHS(runCount),runCount)<FDTHRESH;
                                %Mark one frame before high motion (consistent with Power et al., 2012)
                                oneFrameBefore=find(~threshedMotion)-1;
                                oneFrameBefore=oneFrameBefore(oneFrameBefore>0);
                                threshedMotion(oneFrameBefore)=0;
                                %Mark two frame after high motion (consistent with Power et al., 2012)
                                twoFramesAfter=[find(~threshedMotion)+1 find(~threshedMotion)+2];
                                twoFramesAfter=twoFramesAfter(twoFramesAfter<func_info.RUNLENGTHS(runCount));
                                threshedMotion(twoFramesAfter)=0;
                                %Add high motion points to temporal mask
                                temporalMask(1:func_info.RUNLENGTHS(runCount),runCount)=temporalMask(1:func_info.RUNLENGTHS(runCount),runCount).*threshedMotion;
                                percentTimePointsCensored=100*sum(~threshedMotion)/func_info.RUNLENGTHS(runCount);
                                disp(['Marking high motion time points. ' num2str(percentTimePointsCensored) '% of time points marked for censoring/scrubbing for this run.'])
                                numFramesCensored(runCount)=sum(~threshedMotion);
                            end
                        end
                        runCount=runCount+1;
                    end

                    %**set IMPLEMENT_MOTIONSCRUBBING for REST runs only!
                    if IMPLEMENT_MOTIONSCRUBBING==2
                        nuisanceTSVars.numFramesScrubbedByRun=numFramesCensored;
                        nuisanceTSVars.percentFramesScrubbed=100*sum(numFramesCensored)/func_info.RUNLENGTHS(1);%only applies to Rest
                    end

                    %Organize nuisance regressors
                   nuisanceTSVars.nuisanceTS_wholebrain=nuisanceTS_wholebrain;
                   nuisanceTSVars.temporalMask=temporalMask;
                   
                    nuisanceTSVars.nuisanceTS_whitematter=nuisanceTS_whitematter;
                    nuisanceTSVars.nuisanceTS_ventricles=nuisanceTS_ventricles;
                    nuisanceTSVars.nuisanceTS_motion=nuisanceTS_motion;
                    nuisanceTSVars.FD_motion=FD_motion;
                    nuisanceTSVarsRest=nuisanceTSVars; % use this below for resting state analysis.
                    if modeli>=5 
                        nuisanceTSVars=rmfield(nuisanceTSVars,'nuisanceTS_motion'); %remove 12 regressors of motion
                    end
                   
                    %Save nuisance time series to file (for more efficient processing in future when EXECUTE_PREPNUISANCEREG=0)
                    disp(['Saving results to: ' savedNuisRegfile])
                    save(savedNuisRegfile, 'nuisanceTSVars')
                else
                    %Load nuisance time series from file (for more efficient processing when EXECUTE_PREPNUISANCEREG=0)
                    disp(['Loading results from: ' savedNuisRegfile])
                    % Number of rows = number of TRs from longest run;
                    % Number of cols =number of runs (including rest)
                    load(savedNuisRegfile);        
                end
                output.nuisanceTSVars{subjIndex}=nuisanceTSVars;

                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                %%%% PREPARE TASK REGRESSORS %%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                RERUN_TASKREGRESSORPREP=1;  %Set to 1 if you want to run this procedure again
                savedTaskRegfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_TaskRegressorVars.mat'];
                timingFileThisSub=[timingfileDir 'Viscomp_' subjName '_TaskRegressorVars.mat'];%Only execute this if EXECUTE_TASKREGRESSORPREP==1 or the saved file doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
                if RERUN_TASKREGRESSORPREP==1 % load timing, add the hrf values as a separate structure field, and then resave.
                    if isempty(func_info.TASKRUNS)% || ~exist( timingFileThisSub, 'file')
                       error('Either not enough task runs or the stimulus timing file cannot be found.') 
                    end
                    disp(['Loading results from: '  timingFileThisSub])
                    load(timingFileThisSub); 

                    %Convolve with canonical HRF
                    hrf=spm_hrf(TR_INSECONDS);
                    if modeli==1
                        taskdesignmat=TwoRegresMatBlock;
                    elseif modeli==2
                        taskdesignmat=TwoRegresMat;
                    elseif modeli==3
                        taskdesignmat=FourRegresMatStim;
                    elseif modeli==4
                        taskdesignmat=FourRegresMatRT;
                    elseif modeli==5
                        taskdesignmat=SixRegresMatStimRT; %separate regressors for stim and RT; fewer nuissance regressors
                    elseif modeli==6
                        taskdesignmat=FourRegresMatStim;      % same as Model3, except with fewer nuisance regressors
                    end
                    taskdesignmat_hrf=zeros(size(taskdesignmat));
                    for regressorNum=1:size(taskdesignmat,2)
                        convData=conv(taskdesignmat(:,regressorNum),hrf);
                        taskdesignmat_hrf(:,regressorNum)=convData(1:size(taskdesignmat,1),:);
                    end
                    % define the structure
                    taskTiming.taskdesignmat=taskdesignmat;
                    taskTiming.taskdesignmat_hrf=taskdesignmat_hrf;

                    %Save task timing variables to task regressor file 
                    disp(['Saving results to: ' savedTaskRegfile])
                    save(savedTaskRegfile, 'taskTiming')
                else % file already available
                    if ~isempty(func_info.TASKRUNS)
                        %Load task timing variables from file (for more efficient processing when EXECUTE=0)
                        disp(['Loading results from: ' savedTaskRegfile])
                        load(savedTaskRegfile);
                    end
                end
                if ~isempty(func_info.TASKRUNS)
                    output.taskTiming{subjIndex}=taskTiming;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% REST NUISANCE REGRESSION %%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                RERUN_RESTGLM=1; %Set to 1 if you want to run this procedure again
                savedRestNuisRegfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_RestNuisanceGLMVars.mat'];

                %Only execute this if RERUN_TASKGLM==1 or the savedRestNuisRegfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
                if and(~isempty(func_info.RESTRUNS), or(RERUN_RESTGLM == 1, ~exist(savedRestNuisRegfile, 'file')))
                    disp('Running rest nuisance regression')

                    %Specify the number of nuisance regressors
                    %NUMREGRESSORS_NUISANCE=16;%2 white matter (1 normal + 1 deriv, computed in GLM function below), 2 ventricles (1 normal + 1 deriv, computed in GLM), 12 motion
                    %% TAKU EDIT - Ciric model additions for 36p model (without GSR it's 32; with GSR it's 36 (4 additional regressors for GSR)
                    % 4 white matter (1 normal + 1 derivative, 1 quadratic + 1 derivatives, computed in GLM function below)
                    % 4 ventricles (1 normal + 1 derivative, 1 quadratic + 1 derivative, computed in GLM)
                    % 24 motion (6 normal + 6 derivatives, 6 quadratic + 6 derivatives
                    NUMREGRESSORS_NUISANCE=32;

                    %Add 4 regressors for GSR: 1 normal + 1 derivative, 1 quadratic + 1 derivative
                    if GSR
                        NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE+4;
                    end
                    visualizeDesignMatrix=1;
                    restGLMVars = runGLM(tseriesMatSubj, NUMREGRESSORS_NUISANCE, nuisanceTSVarsRest, [], func_info.RUNLENGTHS, func_info.RESTRUNS, NPROC, NUMPARCELS, GSR, visualizeDesignMatrix);

                   %Save task GLM variables to file (for more efficient processing in future when EXECUTE=0)
                    disp(['Saving results to: ' savedRestNuisRegfile])
                    save(savedRestNuisRegfile, 'restGLMVars')
                else
                    if ~isempty(func_info.RESTRUNS)
                        %Load task GLM variables from file (for more efficient processing when EXECUTE=0)
                        disp(['Loading results from: ' savedRestNuisRegfile])
                        load(savedRestNuisRegfile);
                    end
                end
                if ~isempty(func_info.RESTRUNS)
                    output.rest_fMRI_preprocTS{subjIndex} = restGLMVars.fMRI_resids;
                    output.restGLMVars{subjIndex} = restGLMVars;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% TASK NUISANCE REGRESSION AND GLM %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
                RERUN_TASKGLM=1;
                savedTaskNuisRegGLMfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_TaskGLM.mat'];

                %Only execute this if RERUN_TASKGLM==1 or the savedTaskNuisRegGLMfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
                if and(~isempty(func_info.TASKRUNS), or(RERUN_TASKGLM == 1, ~exist(savedTaskNuisRegGLMfile, 'file')))
                    disp('Running task nuisance regression and GLM')

                    % Specify the number of nuisance regressors
                   NUMREGRESSORS_NUISANCE=32;
                   if modeli>=5
                        NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE-24 %removing 24 motion regressors
                    end

                    % Add 2 regressors for GSR
                    if GSR
                        NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE+4;
                    end
                    visualizeDesignMatrix=1;

                    taskGLMVars = runGLM(tseriesMatSubj, NUMREGRESSORS_NUISANCE, nuisanceTSVars, taskTiming.taskdesignmat_hrf, func_info.RUNLENGTHS, func_info.TASKRUNS, NPROC, NUMPARCELS, GSR, visualizeDesignMatrix);

                    %Save task GLM variables to file (for more efficient processing in future when EXECUTE=0)
                    disp(['Saving results to: ' savedTaskNuisRegGLMfile])
                    save(savedTaskNuisRegGLMfile, 'taskGLMVars')
                else
                    if ~isempty(func_info.TASKRUNS)
                        %Load task GLM variables from file (for more efficient processing when EXECUTE=0)
                        disp(['Loading results from: ' savedTaskNuisRegGLMfile])
                        load(savedTaskNuisRegGLMfile);
                    end
                end

                if ~isempty(func_info.TASKRUNS)
                    output.task_fMRI_preprocTS{subjIndex} = taskGLMVars.fMRI_resids;
                    output.task_betas{subjIndex} = taskGLMVars.fMRI_betas;
                    output.taskGLMVars{subjIndex} = taskGLMVars;
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% TEMPORAL FILTERING (and motion scrubbing) %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
                RERUN_TEMPORALFILTER=1;
                savedTempFiltGLMfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_TemporalFilter.mat'];
                if IMPLEMENT_MOTIONSCRUBBING==2
                    %Only execute this if RERUN_TEMPORALFILTER==1 or the savedTempFiltGLMfile doesn't exist for this subject, but skip if TEMPORALFILTER==0
                    if or(RERUN_TEMPORALFILTER == 1, ~exist(savedTempFiltGLMfile, 'file'))
                        %Create temporal filter (based on Jonathan Power's script from Petersen lab)
                        if TEMPORALFILTER==1
                            disp('Applying temporal filter')
                            lopasscutoff=LOWPASS_HZ/(0.5/TR_INSECONDS); % since TRs vary have to recalc each time
                            hipasscutoff=HIGHPASS_HZ/(0.5/TR_INSECONDS); % since TRs vary have to recalc each time
                            %Using filter order of 1
                            filtorder=1;
                            [butta, buttb]=butter(filtorder,[hipasscutoff lopasscutoff]);
                        end

                        %**set IMPLEMENT_MOTIONSCRUBBING for REST runs only!
                        %Rest data
                        if ~isempty(func_info.RESTRUNS)
                            %Interpolate data if scrubbing
                            %Use interpolation to account for gaps in time series due to motion scrubbing
                            %if IMPLEMENT_MOTIONSCRUBBING==1
                            restfMRIData_scrubbed=restGLMVars.fMRI_resids;
                            fMRIData_rest = interpolateAcrossTSGaps(restfMRIData_scrubbed, restGLMVars.temporalmask, func_info.RUNLENGTHS, func_info.RESTRUNS, NUMPARCELS);
                            %else
                            %    fMRIData_rest=restGLMVars.fMRI_resids;
                            %end

                            %Apply temporal filter
                            if TEMPORALFILTER==1
                                filteredData=filtfilt(butta,buttb,fMRIData_rest');
                                filteredData=filteredData';

                                %Reapply scrubbing
                                if IMPLEMENT_MOTIONSCRUBBING==1
                                    filteredData=filteredData(:,logical(restGLMVars.temporalmask));
                                end
                                filteredDataOutput.rest_fMRI_preprocTS=filteredData;
                            else
                                scrubbedDataOutput.rest_fMRI_preprocTS=fMRIData_rest;
                                %Save filtered time series variables to file (for more efficient processing in future when EXECUTE=0)
                                disp(['Saving results to: ' savedTempFiltGLMfile])
                                save(savedTempFiltGLMfile, 'scrubbedDataOutput')
                            end
                        end
                        %**Not applying filtering to the task runs!
%                             %Task data
%                             if ~isempty(TASKRUNS)
% 
%                                 %Interpolate data if scrubbing
%                                 %Use interpolation to account for gaps in time series due to motion scrubbing
%                                 if IMPLEMENT_MOTIONSCRUBBING==1
%                                     taskfMRIData_scrubbed=taskGLMVars.fMRI_resids;
%                                     fMRIData_task = interpolateAcrossTSGaps(taskfMRIData_scrubbed, taskGLMVars.temporalmask, RUNLENGTHS, TASKRUNS, NUMPARCELS);
%                                 else
%                                     fMRIData_task=taskGLMVars.fMRI_resids;
%                                 end
%
%                                 %Apply temporal filter
%                                 filteredData=filtfilt(butta,buttb,fMRIData_task');
%                                 filteredData=filteredData';
%                                 %Reapply scrubbing
%                                 if IMPLEMENT_MOTIONSCRUBBING==1
%                                     filteredData=filteredData(:,logical(taskGLMVars.temporalmask));
%                                 end
%                                 filteredDataOutput.task_fMRI_preprocTS=filteredData;
%                             end
                    else
                        %Load filtered time series from file (for more efficient processing when EXECUTE=0)
                        disp(['Loading results from: ' savedTempFiltGLMfile])
                        load(savedTempFiltGLMfile);
                    end
%                         if ~isempty(TASKRUNS)
%                             output.task_fMRI_preprocTS{subjIndex}=filteredDataOutput.task_fMRI_preprocTS;
%                         end
                    if ~isempty(func_info.RESTRUNS)
                        output.rest_fMRI_preprocTS_scrubbed{subjIndex}=scrubbedDataOutput.rest_fMRI_preprocTS;
                    end
                end
                close all;
            end
        %Store output in GLMoutput
        GLMOutput.(model_name).(MSM_name)=output;
        output.SUBJECTLIST
       % end  
        end
    msg=sprintf('The model number finished is %d.',modeli);
    warning(msg) % not really a warning, but important to change the font color
    end
end
    
%% Save final result output           
outputfilename=[outputdatadir,ANALYSISNAME,'.mat'];
disp(['Saving final results to: ' outputfilename])
save(outputfilename, 'GLMOutput', '-v7.3');

disp('==Be sure to delete intermediate temporary files when finished with preprocessing all subjects'' data==')
disp('Shell command to delete intermediate temporary files:')
subjTemporaryAnalysisDir_mod=strrep(subjTemporaryAnalysisDir, num2str(SUBJECTLIST{numSubjs}),'*');
disp(['rm -rfv ' subjTemporaryAnalysisDir_mod])
%timeFinished=toc;

%msg2=sprintf('The total time elapsed %d.',modeli) ;
    %disp(msg2)
    
   

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS FOR RUNNING CODE ABOVE %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output_GLM] = runGLM(tseriesMatSubj,NUMREGRESSORS_NUISANCE,nuisanceTSVars,taskdesignmat_hrf,RUNLENGTHS,runNums,NPROC,NUMPARCELS,GSR, visualizeDesignMatrix)
    %Specify the number of task regressors
    numregressors_task=size(taskdesignmat_hrf,2);
    % The two GSR regressors, when called for, count as "NUMREGRESSORS_NUISANCE"
    NUMREGRESSORS=NUMREGRESSORS_NUISANCE+numregressors_task;
    %Add 2 regressors for each run (to account for mean and linear trend within each run)
    numregressors_extra=length(runNums)+length(runNums);

    %Concatenate runs, Organize nuisance regressors
    tseriesMatSubj_fMRIconcat=zeros(NUMPARCELS,sum(RUNLENGTHS(runNums)));
    X=zeros(NUMREGRESSORS+numregressors_extra,sum(RUNLENGTHS(runNums)));
    tmask=ones(sum(RUNLENGTHS(runNums)),1);
    for taskRunIndex=1:length(runNums)
        if taskRunIndex>1
            priorRunsLength=sum(RUNLENGTHS(runNums(1):runNums(taskRunIndex-1)));
        else
            priorRunsLength=0;
        end
        thisRunLength=RUNLENGTHS(runNums(taskRunIndex));
        runStart=priorRunsLength+1;
        runEnd=priorRunsLength+thisRunLength;
        %fMRI data
        tseriesMatSubj_fMRIconcat(:,runStart:runEnd)=tseriesMatSubj(:,1:thisRunLength,runNums(taskRunIndex));
        %regnum=0
        %White matter nuisance regressors
        X(1,runStart:runEnd)=nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex));
        X(2,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex)))];
        % Ciric addition - include quadratics for white matter
        X(3,runStart:runEnd)=nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex)).^2;
        X(4,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex)).^2)];
        %Ventricle
        X(5,runStart:runEnd)=nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex));
        X(6,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)))];
        % Ciric addition - include quadratics for ventricles
        X(7,runStart:runEnd)=nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)).^2;
        %www1=size([0; diff(nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)).^2)])
        %www2=size(diff(nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)).^2))
        X(8,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)).^2)];
        if isfield(nuisanceTSVars,'nuisanceTS_motion')%add motion
            %Motion (12 regressors)
            X(9:20,runStart:runEnd)=nuisanceTSVars.nuisanceTS_motion(:,1:thisRunLength,runNums(taskRunIndex));
            %Quadratic motion regressors (+12 regressors)
            X(21:26,runStart:runEnd)=nuisanceTSVars.nuisanceTS_motion(1:6,1:thisRunLength,runNums(taskRunIndex)).^2;
            %xxx=(size(X(27:32,runStart:runEnd))) % 6 x 765 matrix
            %uuu=size(nuisanceTSVars.nuisanceTS_motion)
            %yyy=size(diff(nuisanceTSVars.nuisanceTS_motion(1:6,1:thisRunLength,runNums(taskRunIndex)), 1, 2).^2) 
            %X(27:32,runStart:runEnd)=[zeros(1,thisRunLength); diff(nuisanceTSVars.nuisanceTS_motion(1:6,1:thisRunLength,runNums(taskRunIndex)).^2)]; % first compute the quadratic, then compute derivative of quadratic
            X(27:32,runStart:runEnd)=[zeros(6,1), diff(nuisanceTSVars.nuisanceTS_motion(1:6,1:thisRunLength,runNums(taskRunIndex)), 1, 2).^2]; 
            %Run global signal regression if specified
            if GSR
                X(33,runStart:runEnd)=nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex));
                X(34,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)))];
                % GSR quadratics
                X(35,runStart:runEnd)=nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)).^2;
                X(36,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)).^2)];
            end
        else
            if GSR
                X(9,runStart:runEnd)=nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex));
                X(10,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)))];
                % GSR quadratics
                X(11,runStart:runEnd)=nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)).^2;
                X(12,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)).^2)];
            end
        end
        
        %Add task regressors
        if numregressors_task>0
            X((NUMREGRESSORS_NUISANCE+1):NUMREGRESSORS,runStart:runEnd)=taskdesignmat_hrf(runStart:runEnd,:)';
        end
        %Run transition regressor
        X(NUMREGRESSORS+taskRunIndex,runStart:runEnd)=ones(thisRunLength,1);
        %Linear trend for run regressor
        X(NUMREGRESSORS+taskRunIndex+length(runNums),runStart:runEnd)=linspace(0,1,thisRunLength);
        %Temporal mask
        tmask(runStart:runEnd)=nuisanceTSVars.temporalMask(1:thisRunLength,runNums(taskRunIndex));
    end

    %Zscore the design matrix to make it easier to visualize
    if visualizeDesignMatrix
        disp('==Make sure to check over the design matrix visually==')
        Xzscored=zscore(X,0,2);
        Xzscored(logical(eye(size(Xzscored))))=0;
        figure;imagesc(Xzscored);title('Regressors');
        disp('Also showing rank correlation among regressors')
        rankMat=zeros(size(X));
        for ind=1:size(X,1)
            rankMat(ind,:)=tiedrank(X(ind,:));
        end
        rankCorrMat=corrcoef(rankMat');
        rankCorrMat(logical(eye(size(rankCorrMat))))=0;
        figure;imagesc(rankCorrMat);title('Regressor Spearman correlations');
    end

    %Apply temporal mask
    X_orig=X;
    X_tmasked=X(:,logical(tmask));
    X=X_tmasked;
    tseriesMatSubj_fMRIconcat=tseriesMatSubj_fMRIconcat(:,logical(tmask));

    %Test rank of design matrix
    matrixRank=rank(X);
    disp(['Number of regressors: ' num2str(size(X,1))])
    disp(['Rank of matrix: ' num2str(matrixRank)])
    if matrixRank < size(X,1)
        disp('ERROR: Matrix is rank deficient; fix the design matrix before running your regression.')
        disp('Consider using PCA on the nuisance regressors (to orthogonalize them). (Do not orthogonalize the task regressors)')
    end

    % Instantiate empty arrays
    fMRI_resids = zeros(size(tseriesMatSubj_fMRIconcat));
    fMRI_betas = zeros(NUMPARCELS, size(X,1));

    % Begin for loop
    parpool(NPROC);
    parfor (regionNum=1:NUMPARCELS, NPROC)
        % Get the region's data
        ROITimeseries = tseriesMatSubj_fMRIconcat(regionNum,:);

        % Regress out the nuisance time series, keep the residuals and betas
        %stats = regstats(ROITimeseries, X', 'linear', {'r', 'beta'});
        [beta,bint,resid] = regress(ROITimeseries', X');

        % Collect rest regression results
        fMRI_resids(regionNum, :) = resid;
        fMRI_betas(regionNum,:) = beta';
    end
    delete(gcp);
    output_GLM.fMRI_resids = fMRI_resids;
    output_GLM.fMRI_betas = fMRI_betas;
    output_GLM.temporalmask = tmask;
    output_GLM.X_orig = X_orig;
    output_GLM.X_tmasked = X_tmasked;
    %Use this for task functional connectivity analyses (e.g., choose time points using a threshold of 0.5)
    if numregressors_task>0
        output_GLM.taskdesignmat_hrf_tmasked = taskdesignmat_hrf(logical(tmask),:);
    end
end


%%%%%%%%%%%%%%%%%%%
function [output_InterpolatedData] = interpolateAcrossTSGaps(fMRI_tseries, tmask, RUNLENGTHS, runNums, NUMPARCELS)
    %Place fMRI data into original-sized matrix, with NaNs in the gaps
    fMRIData_scrubbed=fMRI_tseries;
    fMRIData_withNans=nan(NUMPARCELS,sum(RUNLENGTHS(runNums)));
    fMRIData_withNans(:,logical(tmask))=fMRIData_scrubbed;
    fMRIData_interpolated=zeros(NUMPARCELS,sum(RUNLENGTHS(runNums)));
    for taskRunIndex = 1:length(runNums)
        %Prep run timing and data
        if taskRunIndex>1
            priorRunsLength=sum(RUNLENGTHS(runNums(1):runNums(taskRunIndex-1)));
        else
            priorRunsLength=0;
        end
        thisRunLength=RUNLENGTHS(runNums(taskRunIndex));
        runStart=priorRunsLength+1;
        runEnd=priorRunsLength+thisRunLength;
        %fMRI data
        fMRIdata_thisrun=fMRIData_withNans(:,runStart:runEnd);

        %Identify gaps in time series due to motion scrubbing
        bd=isnan(fMRIdata_thisrun);
        nongap_timepoints=find(~bd);
        bd([1:(min(nongap_timepoints)-1) (max(nongap_timepoints)+1):end])=0;
        %Implement linear interpolation across gaps (but not at beginning and end of run)
        fMRIdata_thisrun_interpolated=fMRIdata_thisrun;
        fMRIdata_thisrun_interpolated(bd)=interp1(nongap_timepoints,fMRIdata_thisrun(nongap_timepoints),find(bd));

        %Set NaNs to 0
        fMRIdata_thisrun_interpolated(isnan(fMRIdata_thisrun_interpolated))=0;

        %Output result
        fMRIData_interpolated(:,runStart:runEnd)=fMRIdata_thisrun_interpolated; 
    end
    output_InterpolatedData=fMRIData_interpolated;
 end



