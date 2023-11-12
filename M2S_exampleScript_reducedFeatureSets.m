%% Test M2S with reduction of number of points of FeatureSets
% First bin with equal bin size, then use only the ones with max(FI) within their bin


% M2S toolbox to match LCMS metabolomics features of untargeted datasets.
% *** Rui Climaco Pinto ***
% *** Imperial College London, 2021 ***


%% Load reference and target features of sets to match 
% These are two data files containing 3 columns without headers, only values for each feature, in this order: RT, MZ, FI. 
% NOTE: RT = retention time median;  MZ = mass-to-charge (m/z) ratio median; FI = Feature Intensity median
% RT should be in minutes. If it is in seconds, divide it by 60;
% If there is no feature intensity use a column with ones or random values.

%% Load the data

datasetName = "LPOS"; % {"LNEG","LPOS"}

% NOTE: there are two Plasma Lipid examples (ESI- and ESI+), each with two
% datasets that can be matched


if datasetName == "LNEG"
    refFilename = 'test_refFeatures_LNEG_csv.csv';
    targetFilename = 'test_targetFeatures_LNEG_csv.csv';
elseif datasetName == "LPOS"
    refFilename = 'test_refFeatures_LPOS_csv.csv';
    targetFilename = 'test_targetFeatures_LPOS_csv.csv';
end

% Function 'importdata.m' loads data from .csv, .txt, or .xlsx 
% (You may need to change directory to the file location)
[refFeatures] = importdata(refFilename);
[targetFeatures] = importdata(targetFilename);

% Create labels for all features

[refMZRT_str] = M2S_createLabelMZRT('ref',refFeatures(:,2),refFeatures(:,1));
[targetMZRT_str] = M2S_createLabelMZRT('target',targetFeatures(:,2),targetFeatures(:,1));

% Visualize the two feature sets
M2S_figureH(0.8,0.5)
subplot(1,2,1), 
M2S_plotMZRT_featureSet(refFeatures,1,8,1); title('Reference featureSet')
subplot(1,2,2), 
M2S_plotMZRT_featureSet(targetFeatures,1,8,1); title('Target featureSet')

%%*************************************************************************
%% Procedure part 1: find all possible matches
%%*************************************************************************

% create a structure to keep the options chosen at each step
opt = struct;

%% Set thresholds for matching all features 
% Try running it with default settings as below, to see the major trends.
% It may become very crowded.
% Default settings: no RTthresh ; MZthresh < 0.02 (m/z); log10FIthresh < 1000
% Later each setting can be set as desired.
[refSet,targetSet,Xr_connIdx,Xt_connIdx, opt]=M2S_matchAll(refFeatures,targetFeatures);
        
% Additional plot of results:
Xsymbol = '.k';
plotRow = 2; % y-axis is e.g. RTdist (1) or both RTdist and RT_target (2)
M2S_plotDelta_matchedSets(refSet,targetSet,Xsymbol,plotRow)


% NOTES on setting thresholds:
% Thresholds define the maximum acceptable inter-dataset shift for each dimension (RT, MZ, log10FI). 
% You can define "horizontal line" fixed thresholds by using only "intercept" values, and setting "slope" to zero.
% These are equivalent to defining a fixed m/z.
% You can define "diagonal line" relative thresholds by using also "slope" values.
% These are equivalent to defining a relative m/z (kind of a ppm threshold value)
% HINT: the function below helps to calculate an intercept and slope in one of the dimensions:
%[RTMZFI_intercept,RTMZFI_slope] = M2S_calculateInterceptSlope() % to define it with the mouse on a figure
%[RTMZFI_intercept,RTMZFI_slope] = M2S_calculateInterceptSlope(startPoint,endPoint) % to insert two points manually

%% M2S_matchAll SETTINGS for the two main examples:

if datasetName == "LNEG"
    % Example settings for LNEG test data
    opt.FIadjustMethod = 'regression'; % {'none','median','regression'}
    opt.multThresh.RT_intercept = [-0.55,0.15];
    opt.multThresh.RT_slope = [0 0];
    opt.multThresh.MZ_intercept = [-0.01 0.01]; % m/z units
    opt.multThresh.MZ_slope = [-5e-6 5e-6]; % ppm
    opt.multThresh.log10FI_intercept = [-1 1.5];
    opt.multThresh.log10FI_slope = [0 1];
    
elseif datasetName == "LPOS"
    % Example settings for LPOS test data
    opt.FIadjustMethod = 'median'; % {'median','regression'}
    opt.multThresh.RT_intercept = [-0.1,0.1];
    opt.multThresh.RT_slope = [0 0];
    opt.multThresh.MZ_intercept = [-0.0075 0.015];
    opt.multThresh.MZ_slope = [-0.01/2000, -0.01/2000];
    opt.multThresh.log10FI_intercept = [-1.4 0.75];
    opt.multThresh.log10FI_slope = [0 0];    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create and work with reduced sets 

% Reduce the size of Feature Sets
[refFeatures_reduced,targetFeatures_reduced, refFeaturesBin_maxFI_idx, targetFeaturesBin_maxFI_idx] = M2S_reduceSet(refFeatures, targetFeatures)

%% Find all matches for the reduced sets

plotType = 3
% [refSetReduced,targetSetReduced,Xr_connIdxReduced,Xt_connIdxReduced, optReduced]=M2S_matchAll(refFeatures_reduced,targetFeatures_reduced);
[refSetReduced,targetSetReduced,Xr_connIdxReduced,Xt_connIdxReduced, optReduced]=M2S_matchAll(refFeatures_reduced,targetFeatures_reduced,opt.multThresh,opt.FIadjustMethod,plotType);

% Calculate the trendline for the reduced dataset
opt.neighbours.nrNeighbors = 0.05;
opt.calculateResiduals.neighMethod = 'cross';
opt.pctPointsLoess = 0.1; % loess with 10% of points
plotTypeResiduals = 1
[Residuals_X_reduced,Reduced_trendline] = M2S_calculateResiduals(refSetReduced,targetSetReduced,Xr_connIdxReduced,Xt_connIdxReduced,opt.neighbours.nrNeighbors, opt.calculateResiduals.neighMethod,opt.pctPointsLoess,plotTypeResiduals)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find all matches for the actual entire sets
plotType = 3; 
[refSet,targetSet,Xr_connIdx,Xt_connIdx,opt]=M2S_matchAll(refFeatures,targetFeatures,opt.multThresh,opt.FIadjustMethod,plotType);

% Find the predicted trendline for all points/matches
[Residuals_X,Residuals_trendline] = M2S_interpolateTrendline(refSetReduced,Reduced_trendline,refSet,targetSet)

% Plot the result to check
figure, 
for p=1:3
subplot(1,3,p)
if p<=2
    plot(refSet(:,p),targetSet(:,p)-refSet(:,p),'.k'), hold on
    plot(refSet(:,p),Residuals_trendline(:,p),'.r')
else
    plot(log10(refSet(:,p)),log10(targetSet(:,p))-log10(refSet(:,p)),'.k'), hold on
    plot(log10(refSet(:,p)),Residuals_trendline(:,p),'.r')
end
axis tight, grid on
end

% Continue with the procedure as usual:
%% Adjust the residuals