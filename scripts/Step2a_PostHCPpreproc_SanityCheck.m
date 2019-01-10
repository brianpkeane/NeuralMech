% Script to ensure that cifti files differ from one another
% Should be 93k brainordinates
% Written by brian p keane in 2018
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Applications/workbench/bin_macosx64/']);
basedir='/Users/briankeane/Desktop/AmarelBackup/data/preprocessed/MNINonLinear/Results/Task_Viscomp';
taskname='Task_Viscomp';
suffix='_Atlas_MSMAll_InitalReg_2_d40_WRN.dtseries.nii';
numRuns=4;
numVertices=91282;
numVols=281;
cii_vc_all=zeros(numVertices,numVols,numRuns);
for i=1:numRuns
    filename=[basedir num2str(i) '/' taskname num2str(i) suffix];
    cii_vc=ciftiopen(filename,'wb_command');
    cii_vc_all(:,:,i)=cii_vc.cdata;
    if i>1
        percentTheSameBetweenAdjacentRuns=100*sum(sum(cii_vc_all(:,:,i)==cii_vc_all(:,:,i-1)))/(numVols*numVertices)
    end
    
end


