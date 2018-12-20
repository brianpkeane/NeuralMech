 
subjName='sub-C02';
BASEDIR='/projects3/NeuralMech/';
datadir=[BASEDIR 'data/preprocessed/'];
timingfileDir=[BASEDIR 'data/stimfiles/']
subjDir=[datadir subjName '/'];
subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
stim_model_input={'m1_2Taskreg','m2_4TaskReg'};
model_name=stim_model_input{1}
MSM_struct_name={'sulcOnly','sulcAll'};
MSM_name=MSM_struct_name{1};
ANALYSISNAME='GLMs_Viscomp';
subjTemporaryAnalysisDir=[subjDir  ANALYSISNAME '/' model_name '_' MSM_name '/']

savedNuisRegfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_nuisanceTSVars.mat'];
savedTaskRegfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_TaskRegressorVars.mat'];
savedRestNuisRegfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_RestNuisanceGLMVars.mat'];
savedTaskNuisRegGLMfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_TaskGLM.mat'];
savedTempFiltGLMfile=[subjTemporaryAnalysisDir 'Viscomp_' subjName '_TemporalFilter.mat'];

