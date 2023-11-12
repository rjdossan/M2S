% Script to match samples and calculate correlation between two datasets
% This is to be used in M2S for the cases when e.g. two different peak
% picking datasets were obtained from the same chromatograms, and we want to match them.
% For example, to match PeakPantheR with xcms output of AIRWAVE cohort.
% Becauset there are actual datasets with the same samples, we can calculate 
% correlation and also see differences in FImed, which are expected to have 
% the exact same values in the two datasets. 
% In the end, instead of the usual penalisation scores, we use the penalisation 
% scores given by the inverted sign correlation values.

% Note: there must be two datasets: Xref_dataset and Xtarget_dataset, each containing (e.g., for Xref_dataset):
% - Xref_dataset.Data : matrix of values only with nRows of samples and nColumns of variables
% - Xref_dataset.SampleInfo : table with nRows (equal to nRows of Data) 
% of sample info, including a variable (column) named "ID", with unique identifiers per sample
% Xref_dataset.VarInfo : table with nRows equal to the nColumns of Data, containing
% info for the variables (metabolites)

%% Create the reference and target datasets

Xref_dataset = XCMS_dataset;
Xref_dataset.SampleInfo.ID = Xref_dataset.SampleInfo.SampleID;% Define

Xtarget_dataset = PPR_dataset;
Xtarget_dataset.SampleInfo.ID = Xtarget_dataset.SampleInfo.SampleID;% Define


%% Match the sample ID of the two datasets

[commonSampleID,sampleID_refIdx,sampleID_targetIdx] = intersect(Xref_dataset.SampleInfo.ID, Xtarget_dataset.SampleInfo.ID);

Xref_dataset.Data = Xref_dataset.Data(sampleID_refIdx,:);
Xref_dataset.SampleInfo = Xref_dataset.SampleInfo(sampleID_refIdx,:);

Xtarget_dataset.Data = Xtarget_dataset.Data(sampleID_targetIdx,:);
Xtarget_dataset.SampleInfo = Xtarget_dataset.SampleInfo(sampleID_targetIdx,:);

%% NOTE: at this point, run the function M2S_matchAll to get Xr_connIdx and Xt_connIdx

% Notice there is FI and the datasets have the same samples, so it is
% possible to calculate correlation between features in the two sets

% Define reference and target feature sets
refFeatures = XCMS_featureSet;
targetFeatures = PPR_featureSet;

% Visualize the two feature sets
M2S_figureH(0.8,0.5)
subplot(1,2,1), 
M2S_plotMZRT_featureSet(refFeatures,1,8,1); title('Reference featureSet')
subplot(1,2,2), 
M2S_plotMZRT_featureSet(targetFeatures,1,8,1); title('Target featureSet')

%%*************************************************************************
%% Procedure part 1: find all possible matches
%%*************************************************************************

% Example of large settings for ROI data

opt.FIadjustMethod = 'median'; % {'median','regression'}
opt.multThresh.RT_intercept = [-0.5,0.5];
opt.multThresh.RT_slope = [0 0];
opt.multThresh.MZ_intercept = [-0.015 0.015];
opt.multThresh.MZ_slope = [0 0];
opt.multThresh.log10FI_intercept = [-100 100];
opt.multThresh.log10FI_slope = [0 0]; 


plotType = 3; 
[refSet,targetSet,Xr_connIdx,Xt_connIdx,opt]=M2S_matchAll(refFeatures,targetFeatures,opt.multThresh,opt.FIadjustMethod,plotType);

%%*************************************************************************
%% Procedure part 2: Calculate penalisation scores for each match
%%*************************************************************************


% Calculate correlation for each of the matches found
corrType = 'Pearson'
corrScore = NaN(length(Xr_connIdx),1);
for matchNr = 1:length(Xr_connIdx)
    corrScore(matchNr) = corr(Xref_dataset.Data(:,Xr_connIdx(matchNr)), Xtarget_dataset.Data(:,Xt_connIdx(matchNr)),'type',corrType);
end

%% Create a penalty score from the corrScore (needs to be inverted)
% Make all negative correlations equal to 1
penaltyScores_corrInverted = -corrScore; % Penalty increases to higher positive values
min(corrScoreFinal)

% [eL,eL_INFO,CC_G1]= M2S_decideBestMatches(refSet, targetSet, Xr_connIdx,Xt_connIdx, penaltyScores); % to compare with correlation score
[eL,eL_INFO,CC_G1]= M2S_decideBestMatches(refSet, targetSet, Xr_connIdx,Xt_connIdx, penaltyScores_corrInverted);

%% Plot the differences in each dimension, coloured by penaltyScores_corrInverted

figure; M2S_plotDelta_matchedSets(refSet,targetSet,'.k',2)
subplot(2,3,1); scatter(refSet(:,1),targetSet(:,1) - refSet(:,1), 20, penaltyScores_corrInverted,'filled'); 
M2S_colorByY_ofSubplot(penaltyScores_corrInverted)
subplot(2,3,1); colorbar
colormap(jet)

% From here continue as usual to obtain the results


