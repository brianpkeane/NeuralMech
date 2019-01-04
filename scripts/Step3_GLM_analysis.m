%Step3_GLM_analysis.m

%OVERALL PLAN (*Focus on Model1a-c comparisons primarily for now; run Model2 if there's time):
% 1. Use ciftiopen to open a dscalar.nii as a template structure: cii = ciftiopen('path/to/file','path/to/wb_command'); CIFTIdata = cii.cdata;
% **Template .dscalar stored as 'EXAMPLE_' in scripts dir; taken from SRActFlow
% 2. Run Mike's GLM script (see NOTES immediately above for parm choices):
% 2a. Parcellate the dtseries file output from HCP preproc using the Glasser parcels (stored in scripts dir). This is done separately for L/R surface hemispheres, and the result concatenated (1-180=L; 181-360=R).
% 2b. Run the GLM - *make sure that StimFiles are read in correctly
%   Will need to run this twice for Practice and Test sessions
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

%TO DO
%1. Verify that orientation is matched; 
%EXAMPLE dscalar is based on dlabel /projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.LR.CorticalAreas_dil_Colors.32k_fs_RL.dlabel.nii'
%whereas parcellation in Step2b was done for (L and R):L_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
%   DONE - seems like if GLM parcellation organizes regions as 1-180=left,
%   181-360=right (as in Mike's script), then the above dlabel file (32k_fs_RL) is in the
%   correction orientation to assign regions to grayordinates. It is ok
%   also to use the Conte69 atlas as an underlay, as the header of the
%   dscalar tells workbench which grayordinates correspond to which
%   hemisphere...
%2. *Need to work on Model2 portion of script (adding averaging across
%betas OR subtraction contrast images to sub_betas portion of script; might
%be good to see for sensory + motor rules?)

%% Set paths to output

addpath('/projects/AnalysisTools/')
addpath('/projects/AnalysisTools/gifti-1.6/')
addpath('/projects/AnalysisTools/ReffuncConverter/')

%SUBJECTLIST = {'sub-1','sub-2','sub-3','sub-6','sub-7','sub-8','sub-9','sub-10','sub-11','sub-13','sub-14','sub-16','sub-17'};

%parcel info
L_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
R_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
parcellationName='Glasser2016';
NUMPARCELS=360;

%Directories
BASEDIR='/projects3/NeuralMech/';
datadir=[BASEDIR '/data/preprocessed/'];
outputdatadir=[BASEDIR '/data/results/GLM/'];

%**Set up Model names (should help input StimFiles + name output struct
stim_model_input={'Model1a_1Taskreg_varblock','Model1b_1Taskreg_consepochlength','Model1c_1Taskreg_varepochRT',...
    'Model2a_64Taskreg_varblock','Model2b_64Taskreg_consepochlength','Model2c_64Taskreg_varepochRT'};
%stim_model_mat={'model1a_reg','model1b_reg','model1c_reg','model2a_reg','model2b_reg','model2c_reg'};

%set input suffix which will be chosen based on loop below (between MSMsulc
%and MSMsulc_MSMall
%MSMsulc output=$sub/MNINonLinear/Results/$RUN/$RUN_Atlas.dtseries.nii
%MSMall output=$sub/MNINonLinear/Results/$RUN/$RUN_Atlas_MSMAll_InitalReg_2_d40_WRN.dtseries.nii
MSM_input={'_Atlas.dtseries.nii','_Atlas_MSMAll_InitalReg_2_d40_WRN.dtseries.nii'};
MSM_struct_name={'MSMsulc','MSMsulc_MSMall'};
session_input={'Prac','Test'};

%*determines GLM output filename
%ANALYSISNAME='GLMs_Model1_Model2';
ANALYSISNAME='GLMs_Model1';

load([outputdatadir,ANALYSISNAME,'.mat']); %GLMoutput

%path to example cii file
%based on this /projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.LR.CorticalAreas_dil_Colors.32k_fs_RL.dlabel.nii'\n",
%Template .dscalar stored as 'EXAMPLE_' in scripts dir; taken from SRActFlow
example_scalar=[BASEDIR '/docs/scripts/EXAMPLE_wholeBrainUnivariateContrast_RMID_v_RINDEX.dscalar.nii']; 
glasser_dlabel='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.LR.CorticalAreas_dil_Colors.32k_fs_RL.dlabel.nii';
    
%**set which modules to execute
run_model1_avg=0;
run_model1_stats=1;

%% 1. Use ciftiopen to open a dscalar.nii as a template structure
%cii = ciftiopen('path/to/file','path/to/wb_command'); CIFTIdata = cii.cdata;

%load up example scalar
cii_scalar=ciftiopen(example_scalar,'wb_command');
CIFTIscalar=cii_scalar.cdata; %contains tstat, pvalue and tstat that survives FDR correction; replace tstat with zstat
num_gray=size(CIFTIscalar,1);

%load up Glasser labels - will need this to replace 
cii_label=ciftiopen(glasser_dlabel,'wb_command');
CIFTIlabel=cii_label.cdata;

%**set thresholds
thresh={'raw','0.05','FDR05'};

%% 2. Group analysis
% Compute contrast images for each subject if necessary (for Model2 only)
% Average betas (or contrast images) for each subject.
% Convert the group tstat to zstat (using andybrainblog code): norminv(tcdf(i1,' num2str(dof) '),0,1), where i1=tstat, dof=numsubs-1
% In AnalysisOutput, store activation maps for zstats thresholded/binarized by significance across following thresholds: p < .05, FDR p < .05, raw.
% *Also store the validation result metric:
    % Plot num voxels with significant activation (based on various p thresholds) = 12 GLMs (model1a,1b,1c,2a,2b,2c x MSMall/MSMsulc) x separately for the 2 prac/test sessions.
    % Plot peak zstats = 12 GLMs x separately for the 2 prac/test sessions.
% Output format:
% AnalysisOutput.(model).(MSMversion).(session).(threshold).beta/zstat/num_sig_voxels/peak_zstat

%loop through GLMoutput
GLM_names=fieldnames(GLMOutput);
%init metrics
sig_regions.Prac.names={};
sig_regions.Prac.peakZ=[];
sig_regions.Prac.beta_std=[];
sig_regions.Prac.p05=[];
sig_regions.Prac.FDR=[];
sig_regions.Test=sig_regions.Prac;

for i=1:length(fieldnames(GLMOutput))
    model_name=GLM_names{i};
    
    %zstat image extraction is different for Model1(a-c) vs Model2(a-c)
    if ~isempty(strfind(model_name,'Model1'))
        %can update for Model2 by providing a vector of col codes; seems
        %like Trial regs go 17,19,21,23 for Prac; 17,19...143 for Test
        trial_col=17;%column in GLMOutput.(model_name).(MSM).(session).task_betas{subjNum,1} with trial beta
        enc_col=18;
        
        %start MSM->session loop; compute beta activation maps
        for j=1:length(MSM_struct_name)
            for k=1:length(session_input)
                
                MSM_name=MSM_struct_name{j};
                session_name=session_input{k};
                
                %set outputname for cifti file
                outname=[model_name,'_',MSM_name,'_',session_name];
                
                %**for model2 would first compute the average over all
                %Trialregs for each subject, before storing this result in
                %sub_betas
                
                %store betas for each subj
                num_subs=length(GLMOutput.(model_name).(MSM_name).(session_name).task_betas);
                sub_betas=[]; %stored in cols
                for s=1:num_subs
                    sub_betas=[sub_betas,GLMOutput.(model_name).(MSM_name).(session_name).task_betas{1,s}(:,trial_col)];
                end
                

                %1. RAW AVERAGE BETA MAP
                execute=run_model1_avg;
                if execute==1
                    avg_beta=mean(sub_betas,2);

                    %convert to function, with inputs
                    %avg_beta,CIFTIlabel,outfile
                    outfile=[outputdatadir,'Rawbeta_',outname,'.dscalar.nii'];
                    %loop through glasssr parcels and 'fill in'
                    %grayordinate rows that correspond to Glasser region y 
                    CIFTIscalar_raw=zeros(num_gray,1);%init output
                    for z=1:size(avg_beta,2) %allows looping of zstat etc
                        for y=1:size(avg_beta,1)
                            glasser_gray=find(CIFTIlabel==y);
                            CIFTIscalar_raw(glasser_gray,z)=avg_beta(y,z);
                        end
                    end

                    %write cifti
                    newcii=cii_scalar; 
                    newcii.cdata=CIFTIscalar_raw;
                    ciftisavereset(newcii,outfile,'wb_command'); %Using ciftisavereset as the data matrices will have a different number of maps/columns from what you started with

                end
                
                %2. RUN STATS - ttest p < .05, ttest FDR p < .05
                %conduct ttest with FDR correction
                %output=zstat,zstat_p05,zstat_FDR05
                execute=run_model1_stats;
                if execute==1
                    [h,p,ci,stats]=ttest(sub_betas,0,'Dim',2);
                    %convert tstat to zstat - %col1 in dscalar output
                    zstat=norminv(tcdf(stats.tstat,num_subs-1),0,1);
                    
                    %save thresholded zstat
                    zstat_p05=zstat;
                    zstat_p05(h==0)=0;
                    
                    %perform FDR correction and save zstat
                    FDR=mafdr(p,'BHFDR',true);
                    zstat_FDR05=zstat;
                    zstat_FDR05(FDR>0.05)=0;
                    
                    %save cifti
                    %write cifti
                    outfile=[outputdatadir,'Zstat_',outname,'.dscalar.nii'];
                    zstat_out=[zstat,zstat_p05,zstat_FDR05];
                    
                    CIFTIscalar_zstat=zeros(num_gray,3);%init output
                    for z=1:size(zstat_out,2) %allows looping of zstat etc
                        for y=1:size(zstat_out,1)
                            glasser_gray=find(CIFTIlabel==y);
                            CIFTIscalar_zstat(glasser_gray,z)=zstat_out(y,z);
                        end
                    end
                    newcii=cii_scalar; 
                    newcii.cdata=CIFTIscalar_zstat;
                    ciftisavereset(newcii,outfile,'wb_command'); %Using ciftisavereset as the data matrices will have a different number of maps/columns from what you started with

                    
                    %also compute metrics:
                    %1.Plot num regions with significant activation (based on various p thresholds) = 12 GLMs (model1a,1b,1c,2a,2b,2c x MSMall/MSMsulc) x separately for the 2 prac/test sessions.
                    p05sig=length(find(zstat_p05~=0));
                    FDRsig=length(find(zstat_FDR05~=0));
                    sig_regions.(session_name).names=[sig_regions.(session_name).names;outname];
                    sig_regions.(session_name).p05=[sig_regions.(session_name).p05;p05sig];
                    sig_regions.(session_name).FDR=[sig_regions.(session_name).FDR;FDRsig];
                    
                    %2.Plot peak zstats = 12 GLMs x separately for the 2 prac/test sessions.
                    sig_regions.(session_name).peakZ=[sig_regions.(session_name).peakZ;max(abs(zstat))];
                    %also compute std for each region across subjects
                    std_beta=std(sub_betas,0,2);
                    sig_regions.(session_name).beta_std=[sig_regions.(session_name).beta_std;mean(std_beta)]; 
                end        
            end
        end
        
    elseif ~isempty(strfind(model_name,'Model2'))
        %compute contrast image for each subject
        
        %average contrast images
        
    end
    
end

%save sig_regions
if run_model1_stats==1
    save([outputdatadir,'ActivationStats.mat'],'sig_regions');
end

% 5. Use workbench to visualize - use the Conte32k atlas as underlay, and the zstat images as the overlay.



