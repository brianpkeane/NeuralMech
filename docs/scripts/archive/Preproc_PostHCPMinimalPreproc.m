function [output] = Preproc_PostHCPMinimalPreproc(SUBJECTLIST)

% This function runs preprocessing and GLM analysis, after the HCP minimal preproccessing pipeline has been run.
% It is designed to parcellate a dataset into a set of regions, which are then preprocessed.
%
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
%
%Script version: 1.5
%Date: October 25, 2016
%
%Version history:
%1.1: Fixed a bug in which temporal filtering would run even if flagged not
%   to. Also changed standard TR duration.
%1.2: Added check for GLM design matrix deficiency
%1.3: NPROC now works properly. Added information on deleting intermediate 
%   temporary files to save disk space. Changed which directory files are 
%   saved to (now analysis_[ANALYSISNAME]). Now saves out 
%   output_GLM.taskdesignmat_hrf_tmasked, which indicates task time points 
%   after scrubbing.
%1.4: Fixed bug in which temporal filtering was only run on the last 
%   subject (instead of all subjects). Changed final output file to v7.3 to
%   accomodate large file sizes.
%1.5: Fixed a bug in which the hemispheres might have been flipped. Now each
%   hemisphere is loaded separately to ensure they are not flipped.
%1.6: Added new paths to the script to include more dependencies

%% Parameters to customize for your analysis
addpath('/projects/AnalysisTools/')
addpath('/projects/AnalysisTools/gifti-1.6/')
addpath('/projects/AnalysisTools/ReffuncConverter/')
%Dataset-specific paths
addpath('/projects3/NeuralMech/docs/scripts/')

%%Basic processing parameters
ANALYSISNAME='TaskFCMethods1';
%subjList = list of subject numbers as a list of strings (used if not set by function call)
if isempty(SUBJECTLIST)
    SUBJECTLIST = {'sub-S31'};
    %Full subject list: {'013', '014', '016', '017', '018', '021', '023', '024', '025', '026', '027', '028', '030', '031', '032', '033', '034', '035', '037', '038', '039', '040', '041', '042', '043', '045', '046', '047', '048', '049', '050', '053', '055', '056', '057', '058', '062', '063', '064', '066', '067', '068', '069', '070', '072', '074', '075', '076', '077', '081', '082', '085', '086', '087', '088', '090', '092', '093', '094', '095', '097', '098', '099', '101', '102', '103', '104', '105', '106', '108', '109', '110', '111', '112', '114', '115', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '134', '135', '136', '137', '138', '139', '140', '141'};
end
numSubjs=length(SUBJECTLIST);
TR_INSECONDS=0.785;
FRAMESTOSKIP=5;

%Basic data parameters
RUNNAMES = {'Task_Rest','Task_Retino1','Task_Retino2','Task_Retino3','Task_Viscomp1','Task_Viscomp2','Task_Viscomp3','Task_Viscomp4'};
numRuns=length(RUNNAMES);
RESTRUNS=1;
TASKRUNS=2:8;
RUNLENGTHS = [765, 306, 306, 306, 281, 281, 281, 281];
L_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
R_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
parcellationName='Glasser2016';
NUMPARCELS=360;

%Data processing flags
GSR=0;      %GSR = 0 if no GSR, 1 if you want to include GSR
NPROC=4;    %Number of processors to use when computing GLMs
IMPLEMENT_MOTIONSCRUBBING=1;     %Set to 1 for yes, 0 for no
FDTHRESH=0.25;  %Framewise displacement; The threshold (in millimeters) for flagging in-scanner movement
TEMPORALFILTER=1;   %Flag indicating if data should be temporally filtered
HIGHPASS_HZ=0.008;
LOWPASS_HZ=0.09;

%Directories
BASEDIR='/projects3/NeuralMech/';
datadir=[BASEDIR '/data/HCPoutput'];
outputdatadir='/projects3/NeuralMech/data/fullypreprocdata/';
timingfileDir=[datadir 'timingfiles/'];
if ~exist(timingfileDir, 'dir'); mkdir(timingfileDir); end

%% Iterate through subjects

output=[];
output.SUBJECTLIST=SUBJECTLIST;
output.RUNNAMES=RUNNAMES;
output.RUNLENGTHS=RUNLENGTHS;
output.RESTRUNS=RESTRUNS;
output.TASKRUNS=TASKRUNS;
output.parcellationName=parcellationName;

for subjIndex=1:numSubjs
    subjNum = SUBJECTLIST{subjIndex};
    disp(['Processing subject ' subjNum]);
    
    tseriesMatSubj=zeros(NUMPARCELS,max(RUNLENGTHS),numRuns);
    
    %Directories
    subjDir=[datadir '/' subjNum '/'];
    SUBJDIROUTPUT=[datadir '/' subjNum '/'];   %Typically set to be same as subjDir
    if ~exist(SUBJDIROUTPUT, 'dir'); mkdir(SUBJDIROUTPUT); end
    subjTemporaryAnalysisDir=[SUBJDIROUTPUT '/analysis_' ANALYSISNAME '/'];
    if ~exist(subjTemporaryAnalysisDir, 'dir'); mkdir(subjTemporaryAnalysisDir); end
    
    
    %% Downsampling grayordinate data to parcels
    
    disp('Downsampling grayordinate data to parcels')
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_EXTRACTPARCELS=0;
    
    runCount=1;
    for runName=RUNNAMES
        thisRunName=runName{1};
        
        L_parcelTSFilename=[subjTemporaryAnalysisDir '/' thisRunName '_Atlas.L.' parcellationName 'Parcels.32k_fs_LR.ptseries.nii'];
        R_parcelTSFilename=[subjTemporaryAnalysisDir '/' thisRunName '_Atlas.R.' parcellationName 'Parcels.32k_fs_LR.ptseries.nii'];
        
        if or(RERUN_EXTRACTPARCELS == 1, ~exist(parcelTSFilename, 'file'))
            
            subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
            
            inputFile=[subjRunDir thisRunName '_Atlas.dtseries.nii'];
            
            eval(['!wb_command -cifti-parcellate ' inputFile ' ' L_parcelCIFTIFile ' COLUMN ' L_parcelTSFilename ' -method MEAN'])
            eval(['!wb_command -cifti-parcellate ' inputFile ' ' R_parcelCIFTIFile ' COLUMN ' R_parcelTSFilename ' -method MEAN'])
            
        end
        
        %Load parcellated data
        L_dat = ciftiopen(L_parcelTSFilename,'wb_command');
        R_dat = ciftiopen(R_parcelTSFilename,'wb_command');
        if size(L_dat.cdata,2)>RUNLENGTHS(runCount)
            disp(['WARNING: More TRs for this run than expected. Subject: ' subjNum ', Run: ' num2str(runCount)])
        end
        tseriesMatSubj(1:180,1:RUNLENGTHS(runCount),runCount)=L_dat.cdata(:,1:RUNLENGTHS(runCount));
        tseriesMatSubj(181:end,1:RUNLENGTHS(runCount),runCount)=R_dat.cdata(:,1:RUNLENGTHS(runCount));
        
        runCount=runCount+1;
    end

    
    %% Prepare nuisance regressors
    disp('Preparing nuisance regressors')
    %Using Freesurfer aparc+aseg masks
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_PREPNUISANCEREG=0;
    
    savedfile=[subjTemporaryAnalysisDir mfilename '_' ANALYSISNAME '_nuisanceTSVars.mat'];
    
    if or(RERUN_PREPNUISANCEREG == 1, ~exist(savedfile, 'file'))
        
        subjMaskDir=[subjDir '/masks/'];
        if ~exist(subjMaskDir, 'dir'); mkdir(subjMaskDir); end
        
        %Resample Freesurfer segmented mask into functional space (using nearest neighbor interpolation); uses AFNI
        subjMNINonLinearDir=['/projects3/NeuralMech/data/HCPoutput/' subjNum '/MNINonLinear/'];
        exampleFunctionalVolFile=[subjDir '/MNINonLinear/Results/' RUNNAMES{1} '/' RUNNAMES{1} '.nii.gz'];
        eval(['!3dresample -overwrite -rmode NN -master ' exampleFunctionalVolFile ' -inset ' subjMNINonLinearDir 'aparc+aseg.nii.gz -prefix ' subjMNINonLinearDir 'aparc+aseg_resampFunc.nii.gz']);
        
        %Load Freesurfer segmented mask
        aparc_aseg=load_nifti(['/projects3/NeuralMech/data/HCPoutput/' subjNum '/MNINonLinear/aparc+aseg_resampFunc.nii.gz']);
        
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
        nuisanceTS_whitematter=zeros(max(RUNLENGTHS),numRuns);
        nuisanceTS_ventricles=zeros(max(RUNLENGTHS),numRuns);
        nuisanceTS_wholebrain=zeros(max(RUNLENGTHS),numRuns);
        nuisanceTS_motion=zeros(12,max(RUNLENGTHS),numRuns);
        FD_motion=zeros(max(RUNLENGTHS),numRuns);
        temporalMask=ones(max(RUNLENGTHS),numRuns);
        numFramesCensored=zeros(numRuns,1);
        
        runCount=1;
        for runName=RUNNAMES
            thisRunName=runName{1};
            
            subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
            
            inputFile=[subjRunDir thisRunName '.nii.gz'];
            
            %Load data
            runData=load_nifti(inputFile);
            runData2D=reshape(runData.vol,size(runData.vol,1)*size(runData.vol,2)*size(runData.vol,3),size(runData.vol,4));
            if size(runData2D,2)>RUNLENGTHS(runCount)
                disp(['WARNING: More TRs for this run than expected. Subject: ' subjNum ', Run: ' thisRunName ', Run number: ' num2str(runCount)])
                runData2D=runData2D(:,1:RUNLENGTHS(runCount));
            end
            
            whiteMatterMask_eroded_1D=reshape(whiteMatterMask_eroded,size(whiteMatterMask_eroded,1)*size(whiteMatterMask_eroded,2)*size(whiteMatterMask_eroded,3),1);
            nuisanceTS_whitematter(1:RUNLENGTHS(runCount),runCount)=mean(runData2D(whiteMatterMask_eroded_1D,:),1);
            
            ventricleMask_eroded_1D=reshape(ventricleMask_eroded,size(ventricleMask_eroded,1)*size(ventricleMask_eroded,2)*size(ventricleMask_eroded,3),1);
            nuisanceTS_ventricles(1:RUNLENGTHS(runCount),runCount)=mean(runData2D(ventricleMask_eroded_1D,:),1);
            
            wholebrainMask_1D=reshape(wholebrainMask,size(wholebrainMask,1)*size(wholebrainMask,2)*size(wholebrainMask,3),1);
            nuisanceTS_wholebrain(1:RUNLENGTHS(runCount),runCount)=mean(runData2D(wholebrainMask_1D,:),1);
            
            %Note: derivatives are already included in motion time series
            motionvals=importdata([subjRunDir 'Movement_Regressors.txt'])';
            nuisanceTS_motion(:,1:RUNLENGTHS(runCount),runCount)=motionvals(:,1:RUNLENGTHS(runCount));
            
            %Skip first FRAMESTOSKIP frames
            temporalMask(1:FRAMESTOSKIP,runCount)=0;
            
            %Calculate framewise displacement (FD) according to Power et al. (2012)
            %Briefly: The sum of the absolute values of the translational and rotational displacements over all frames (in mm)
            %Note: HCP's minimal preprocessing pipeline uses the following ordering of the motion parameters (see https://github.com/Washington-University/Pipelines/blob/master/global/scripts/mcflirt_acc.sh):
            %trans x, trans y, trans z, rot x, rot y, rot z [rotations in degrees], then derivatives of those 6 (for 12 total)
            motionTS_dt=[zeros(size(nuisanceTS_motion,1),1) diff(squeeze(nuisanceTS_motion(:,1:RUNLENGTHS(runCount),runCount))')'];
            assumedRadius=50;
            rot_x=(2*assumedRadius*pi/360)*motionTS_dt(4,:);
            rot_y=(2*assumedRadius*pi/360)*motionTS_dt(5,:);
            rot_z=(2*assumedRadius*pi/360)*motionTS_dt(6,:);
            FD_motion(1:RUNLENGTHS(runCount),runCount)=abs(motionTS_dt(1,:))+abs(motionTS_dt(2,:))+abs(motionTS_dt(3,:))+abs(rot_x)+abs(rot_y)+abs(rot_z);
            %Apply temporal filtering to FD, to reduce the effect of respiration on FD measure. Based on Siegel et al. (2016) [Siegel JS, Mitra A, Laumann TO, Seitzman BA, Raichle M, Corbetta M, Snyder AZ (2016) ?Data Quality Influences Observed Links Between Functional Connectivity and Behavior?. Cereb Cortex. 1?11.http://doi.org/10.1093/cercor/bhw253]
            %Create temporal filter
            lopasscutoff=0.3/(0.5/TR_INSECONDS); % lowpass filter of 0.3 Hz
            %Using filter order of 1
            filtorder=1;
            [butta, buttb]=butter(filtorder,lopasscutoff);
            %Apply temporal filter
            filteredFD=filtfilt(butta,buttb,FD_motion(1:RUNLENGTHS(runCount),runCount));
            FD_motion(1:RUNLENGTHS(runCount),runCount)=filteredFD;
            
            %Implement motion scrubbing/censoring in the temporal mask
            if IMPLEMENT_MOTIONSCRUBBING
                threshedMotion=FD_motion(1:RUNLENGTHS(runCount),runCount)<FDTHRESH;
                %Mark one frame before high motion (consistent with Power et al., 2012)
                oneFrameBefore=find(~threshedMotion)-1;
                oneFrameBefore=oneFrameBefore(oneFrameBefore>0);
                threshedMotion(oneFrameBefore)=0;
                %Mark two frame after high motion (consistent with Power et al., 2012)
                twoFramesAfter=[find(~threshedMotion)+1 find(~threshedMotion)+2];
                twoFramesAfter=twoFramesAfter(twoFramesAfter<RUNLENGTHS(runCount));
                threshedMotion(twoFramesAfter)=0;
                %Add high motion points to temporal mask
                temporalMask(1:RUNLENGTHS(runCount),runCount)=temporalMask(1:RUNLENGTHS(runCount),runCount).*threshedMotion;
                percentTimePointsCensored=100*sum(~threshedMotion)/RUNLENGTHS(runCount);
                disp(['Marking high motion time points. ' num2str(percentTimePointsCensored) '% of time points marked for censoring/scrubbing for this run.'])
                numFramesCensored(runCount)=sum(~threshedMotion);
            end
            
            runCount=runCount+1;
        end
        
        if IMPLEMENT_MOTIONSCRUBBING
            nuisanceTSVars.numFramesScrubbedByRun=numFramesCensored;
            nuisanceTSVars.percentFramesScrubbed=100*sum(numFramesCensored)/sum(RUNLENGTHS);
        end
        
        %Organize nuisance regressors
        nuisanceTSVars.nuisanceTS_whitematter=nuisanceTS_whitematter;
        nuisanceTSVars.nuisanceTS_ventricles=nuisanceTS_ventricles;
        nuisanceTSVars.nuisanceTS_wholebrain=nuisanceTS_wholebrain;
        nuisanceTSVars.nuisanceTS_motion=nuisanceTS_motion;
        nuisanceTSVars.FD_motion=FD_motion;
        nuisanceTSVars.temporalMask=temporalMask;
        
        %Save nuisance time series to file (for more efficient processing in future when EXECUTE_PREPNUISANCEREG=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'nuisanceTSVars')
        
    else
        
        %Load nuisance time series from file (for more efficient processing when EXECUTE_PREPNUISANCEREG=0)
        disp(['Loading results from: ' savedfile])
        load(savedfile);        
        
    end
    
    output.nuisanceTSVars{subjIndex}=nuisanceTSVars;
    
    
    
    %% Prepare task regressors
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_TASKREGRESSORPREP=1;
    
    savedfile=[subjTemporaryAnalysisDir mfilename '_' ANALYSISNAME '_TaskRegressorVars.mat'];
    
    %Only execute this if EXECUTE_TASKREGRESSORPREP==1 or the savedfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
    if and(~isempty(TASKRUNS), or(RERUN_TASKREGRESSORPREP == 1, ~exist(savedfile, 'file')))
        
        %% PROJECT-SPECIFIC CUSTOM SCRIPTS THAT IMPORT BEHAVIORAL DATA INFO INTO MATLAB
        % The end result should be a "reference function"; a list of strings labeling each time point (TR) with behavioral info
        
        % Import fMRI Behavioral EDAT files
        EDATIMPORT = EDATImportBySubj_IndivRITL('_fMRI_CPRO.txt', {subjNum});
        % Create reference function, converting trial information to TR-by-TR info
        % Note: "Good" miniblock means at least 2/3 correct, otherwise marked "bad"
        reffunc = takuReffunc_ModalityControl_v2(EDATIMPORT, 'full64tasks_TaskFCMethods');
        % Put all subjects into a single vector
        %reffunc_vector = takuReffuncVector_IndivRITL(reffunc);
        
        %% Convert reference function to regressors, in AFNI format (can also be used with MATLAB regression functions)
        % ReffuncToAFNI - create stimulus timingfiles in the current directory
        [taskdesignmat, dLabels] = ReffuncToAFNI(reffunc{1}.RefFunc, {subjNum}, 7, [1:7], [306 306 306 281 281 281 281], TR_INSECONDS,...
            {'Task1', ...
            'Task2', ...
            'Task3', ...
            'Task4', ...
            'Task5', ...
            'Task6', ...
            'Task7', ...
            'Task8', ...
            'Task9', ...
            'Task10', ...
            'Task11', ...
            'Task12', ...
            'Task13', ...
            'Task14', ...
            'Task15', ...
            'Task16', ...
            'Task17', ...
            'Task18', ...
            'Task19', ...
            'Task20', ...
            'Task21', ...
            'Task22', ...
            'Task23', ...
            'Task24', ...
            'Task25', ...
            'Task26', ...
            'Task27', ...
            'Task28', ...
            'Task29', ...
            'Task30', ...
            'Task31', ...
            'Task32', ...
            'Task33', ...
            'Task34', ...
            'Task35', ...
            'Task36', ...
            'Task37', ...
            'Task38', ...
            'Task39', ...
            'Task40', ...
            'Task41', ...
            'Task42', ...
            'Task43', ...
            'Task44', ...
            'Task45', ...
            'Task46', ...
            'Task47', ...
            'Task48', ...
            'Task49', ...
            'Task50', ...
            'Task51', ...
            'Task52', ...
            'Task53', ...
            'Task54', ...
            'Task55', ...
            'Task56', ...
            'Task57', ...
            'Task58', ...
            'Task59', ...
            'Task60', ...
            'Task61', ...
            'Task62', ...
            'Task63', ...
            'Task64'},...
            {'\w*_TaskNum1_\w*', ...
            '\w*_TaskNum2_\w*', ...
            '\w*_TaskNum3_\w*', ...
            '\w*_TaskNum4_\w*', ...
            '\w*_TaskNum5_\w*', ...
            '\w*_TaskNum6_\w*', ...
            '\w*_TaskNum7_\w*', ...
            '\w*_TaskNum8_\w*', ...
            '\w*_TaskNum9_\w*', ...
            '\w*_TaskNum10_\w*', ...
            '\w*_TaskNum11_\w*', ...
            '\w*_TaskNum12_\w*', ...
            '\w*_TaskNum13_\w*', ...
            '\w*_TaskNum14_\w*', ...
            '\w*_TaskNum15_\w*', ...
            '\w*_TaskNum16_\w*', ...
            '\w*_TaskNum17_\w*', ...
            '\w*_TaskNum18_\w*', ...
            '\w*_TaskNum19_\w*', ...
            '\w*_TaskNum20_\w*', ...
            '\w*_TaskNum21_\w*', ...
            '\w*_TaskNum22_\w*', ...
            '\w*_TaskNum23_\w*', ...
            '\w*_TaskNum24_\w*', ...
            '\w*_TaskNum25_\w*', ...
            '\w*_TaskNum26_\w*', ...
            '\w*_TaskNum27_\w*', ...
            '\w*_TaskNum28_\w*', ...
            '\w*_TaskNum29_\w*', ...
            '\w*_TaskNum30_\w*', ...
            '\w*_TaskNum31_\w*', ...
            '\w*_TaskNum32_\w*', ...
            '\w*_TaskNum33_\w*', ...
            '\w*_TaskNum34_\w*', ...
            '\w*_TaskNum35_\w*', ...
            '\w*_TaskNum36_\w*', ...
            '\w*_TaskNum37_\w*', ...
            '\w*_TaskNum38_\w*', ...
            '\w*_TaskNum39_\w*', ...
            '\w*_TaskNum40_\w*', ...
            '\w*_TaskNum41_\w*', ...
            '\w*_TaskNum42_\w*', ...
            '\w*_TaskNum43_\w*', ...
            '\w*_TaskNum44_\w*', ...
            '\w*_TaskNum45_\w*', ...
            '\w*_TaskNum46_\w*', ...
            '\w*_TaskNum47_\w*', ...
            '\w*_TaskNum48_\w*', ...
            '\w*_TaskNum49_\w*', ...
            '\w*_TaskNum50_\w*', ...
            '\w*_TaskNum51_\w*', ...
            '\w*_TaskNum52_\w*', ...
            '\w*_TaskNum53_\w*', ...
            '\w*_TaskNum54_\w*', ...
            '\w*_TaskNum55_\w*', ...
            '\w*_TaskNum56_\w*', ...
            '\w*_TaskNum57_\w*', ...
            '\w*_TaskNum58_\w*', ...
            '\w*_TaskNum59_\w*', ...
            '\w*_TaskNum60_\w*', ...
            '\w*_TaskNum61_\w*', ...
            '\w*_TaskNum62_\w*', ...
            '\w*_TaskNum63_\w*', ...
            '\w*_TaskNum64_\w*'},...
            1, 0, 'Unique64TaskSet_IndivRITL');
        
        %Move resulting timing files to timing files directory (if saving to AFNI .1D files)
        %eval(['!mv *stimfile*.1D ' timingfileDir]);
        %eval(['!mv *ConcatRef*.1D ' timingfileDir]);
        
        %Convolve with canonical hemodynamic response function (HRF)
        hrf=spm_hrf(TR_INSECONDS);
        taskdesignmat_hrf=zeros(size(taskdesignmat));
        for regressorNum=1:size(taskdesignmat,2)
            convData=conv(taskdesignmat(:,regressorNum),hrf);
            taskdesignmat_hrf(:,regressorNum)=convData(1:size(taskdesignmat,1),:);
        end
        
        taskTiming.taskdesignmat=taskdesignmat;
        taskTiming.taskdesignmat_hrf=taskdesignmat_hrf;
        taskTiming.EDATIMPORT=EDATIMPORT;
        taskTiming.reffunc=reffunc;
        
        %Save task timing variables to file (for more efficient processing in future when EXECUTE=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'taskTiming')
        
    else
        
        if ~isempty(TASKRUNS)
            %Load task timing variables from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
        end
        
    end
    
    if ~isempty(TASKRUNS)
        output.taskTiming{subjIndex}=taskTiming;
    end
    
    
    %% Rest nuisance regression
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_RESTGLM=1;
    
    savedfile=[subjTemporaryAnalysisDir mfilename '_' ANALYSISNAME '_RestNuisanceGLMVars.mat'];
    
    %Only execute this if RERUN_TASKGLM==1 or the savedfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
    if and(~isempty(RESTRUNS), or(RERUN_RESTGLM == 1, ~exist(savedfile, 'file')))
        disp('Running rest nuisance regression')
        
        %Specify the number of nuisance regressors
        NUMREGRESSORS_NUISANCE=16;
        %Add 2 regressors for GSR
        if GSR
            NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE+2;
        end
        visualizeDesignMatrix=0;
        
        restGLMVars = runGLM(tseriesMatSubj, NUMREGRESSORS_NUISANCE, nuisanceTSVars, [], RUNLENGTHS, RESTRUNS, NPROC, NUMPARCELS, GSR, visualizeDesignMatrix);
        
        
       %Save task GLM variables to file (for more efficient processing in future when EXECUTE=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'restGLMVars')
        
    else
        
        if ~isempty(RESTRUNS)
            %Load task GLM variables from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
        end
        
    end
    
    if ~isempty(RESTRUNS)
        output.rest_fMRI_preprocTS{subjIndex} = restGLMVars.fMRI_resids;
        output.restGLMVars{subjIndex} = restGLMVars;
    end
    
    
    %% Task nuisance regression and GLM
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_TASKGLM=1;
    
    savedfile=[subjTemporaryAnalysisDir mfilename '_' ANALYSISNAME '_TaskGLM.mat'];
    
    %Only execute this if RERUN_TASKGLM==1 or the savedfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
    if and(~isempty(TASKRUNS), or(RERUN_TASKGLM == 1, ~exist(savedfile, 'file')))

        disp('Running task nuisance regression and GLM')
        
        %Specify the number of nuisance regressors
        NUMREGRESSORS_NUISANCE=16;
        %Add 2 regressors for GSR
        if GSR
            NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE+2;
        end
        visualizeDesignMatrix=0;
        
        taskGLMVars = runGLM(tseriesMatSubj, NUMREGRESSORS_NUISANCE, nuisanceTSVars, taskTiming.taskdesignmat_hrf, RUNLENGTHS, TASKRUNS, NPROC, NUMPARCELS, GSR, visualizeDesignMatrix);
        
        %Save task GLM variables to file (for more efficient processing in future when EXECUTE=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'taskGLMVars')
        
    else
        
        if ~isempty(TASKRUNS)
            %Load task GLM variables from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
        end
        
    end
    
    if ~isempty(TASKRUNS)
        output.task_fMRI_preprocTS{subjIndex} = taskGLMVars.fMRI_resids;
        output.task_betas{subjIndex} = taskGLMVars.fMRI_betas;
        output.taskGLMVars{subjIndex} = taskGLMVars;
    end
    
    
    %% Apply temporal filtering
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_TEMPORALFILTER=1;
    
    savedfile=[subjTemporaryAnalysisDir mfilename '_' ANALYSISNAME '_TemporalFilter.mat'];
    
    if TEMPORALFILTER==1
        
        %Only execute this if RERUN_TEMPORALFILTER==1 or the savedfile doesn't exist for this subject, but skip if TEMPORALFILTER==0
        if or(RERUN_TEMPORALFILTER == 1, ~exist(savedfile, 'file'))
            
            disp('Applying temporal filter')
            
            %Create temporal filter (based on Jonathan Power's script from Petersen lab)
            lopasscutoff=LOWPASS_HZ/(0.5/TR_INSECONDS); % since TRs vary have to recalc each time
            hipasscutoff=HIGHPASS_HZ/(0.5/TR_INSECONDS); % since TRs vary have to recalc each time
            %Using filter order of 1
            filtorder=1;
            [butta, buttb]=butter(filtorder,[hipasscutoff lopasscutoff]);
            
            %Rest data
            if ~isempty(RESTRUNS)
                
                %Interpolate data if scrubbing
                %Use interpolation to account for gaps in time series due to motion scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    restfMRIData_scrubbed=restGLMVars.fMRI_resids;
                    fMRIData_rest = interpolateAcrossTSGaps(restfMRIData_scrubbed, restGLMVars.temporalmask, RUNLENGTHS, RESTRUNS, NUMPARCELS);
                else
                    fMRIData_rest=restGLMVars.fMRI_resids;
                end
                
                %Apply temporal filter
                filteredData=filtfilt(butta,buttb,fMRIData_rest');
                filteredData=filteredData';
                
                %Reapply scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    filteredData=filteredData(:,logical(restGLMVars.temporalmask));
                end
                
                filteredDataOutput.rest_fMRI_preprocTS=filteredData;
                
            end
            
            %Task data
            if ~isempty(TASKRUNS)
                
                %Interpolate data if scrubbing
                %Use interpolation to account for gaps in time series due to motion scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    taskfMRIData_scrubbed=taskGLMVars.fMRI_resids;
                    fMRIData_task = interpolateAcrossTSGaps(taskfMRIData_scrubbed, taskGLMVars.temporalmask, RUNLENGTHS, TASKRUNS, NUMPARCELS);
                else
                    fMRIData_task=taskGLMVars.fMRI_resids;
                end
                
                %Apply temporal filter
                filteredData=filtfilt(butta,buttb,fMRIData_task');
                filteredData=filteredData';
                
                %Reapply scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    filteredData=filteredData(:,logical(taskGLMVars.temporalmask));
                end
                
                filteredDataOutput.task_fMRI_preprocTS=filteredData;
                
            end
            
            %Save filtered time series variables to file (for more efficient processing in future when EXECUTE=0)
            disp(['Saving results to: ' savedfile])
            save(savedfile, 'filteredDataOutput')
            
        else
            
            %Load filtered time series from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
            
        end
        
        if ~isempty(TASKRUNS)
            output.task_fMRI_preprocTS{subjIndex}=filteredDataOutput.task_fMRI_preprocTS;
        end
        if ~isempty(RESTRUNS)
            output.rest_fMRI_preprocTS{subjIndex}=filteredDataOutput.rest_fMRI_preprocTS;
        end
        
    end
end


%% Save final result output

outputfilename=[outputdatadir 'results/' mfilename '_' ANALYSISNAME '_output.mat'];
disp(['Saving final results to: ' outputfilename])
save(outputfilename, 'output', '-v7.3');

disp('==Be sure to delete intermediate temporary files when finished with preprocessing all subjects'' data==')
disp('Shell command to delete intermediate temporary files:')
subjTemporaryAnalysisDir_mod=strrep(subjTemporaryAnalysisDir, num2str(SUBJECTLIST{numSubjs}),'*');
disp(['rm -rfv ' subjTemporaryAnalysisDir_mod])

end


%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output_GLM] = runGLM(tseriesMatSubj,NUMREGRESSORS_NUISANCE,nuisanceTSVars,taskdesignmat_hrf,RUNLENGTHS,runNums,NPROC,NUMPARCELS,GSR, visualizeDesignMatrix)

%Specify the number of task regressors
numregressors_task=size(taskdesignmat_hrf,2);
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
    %White matter nuisance regressors
    X(1,runStart:runEnd)=nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex));
    X(2,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex)))];
    %Ventricle
    X(3,runStart:runEnd)=nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex));
    X(4,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)))];
    %Motion (12 regressors)
    X(5:16,runStart:runEnd)=nuisanceTSVars.nuisanceTS_motion(:,1:thisRunLength,runNums(taskRunIndex));
    %Run global signal regression if specified
    if GSR
        X(17,runStart:runEnd)=nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex));
        X(18,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)))];
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



