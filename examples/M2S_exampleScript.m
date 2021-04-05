%% Main script example for matching two untargeted LCMS datasets
% This script contains brief explanations on each step of the M2S toolbox to match metabolomics features of two untargeted datasets using RT/MZ/FI.
% It contains data to run two examples, all possibilities of inputs, and lines of code to run to run the method on default and another on defined settings.
% *** To run code and create detached figures select the desired line(s) and click F9 ***
% You can save this file with another name and then modify it to run your own experiments.

% *** Rui Climaco Pinto ***
% *** Imperial College London, 2021 ***

%% Load reference and target features of sets to match 
% These are files containing 3 columns without headers, only values for each feature, in this order: RT, MZ, FI. 
% NOTE: RT = retention time median;  MZ = mass-to-charge (m/z) ratio median; FI = Feature Intensity median
% RT should be in minutes. If it is in seconds, divide it by 60;
% If there is no feature intensity use a column with ones or random values.


% NOTE: there are two Plasma Lipid examples below, one in NEG and the other in POS mode
% Only ref and target datasets in the same mode can be matched.
refFilename = 'test_refFeatures_LNEG_csv.csv';
% refFilename = 'test_refFeatures_LPOS_csv.csv';
targetFilename = 'test_targetFeatures_LNEG_csv.csv';
% targetFilename = 'test_targetFeatures_LPOS_csv.csv';

% Function 'importdata.m' loads data from .csv, .txt, or .xlsx 

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
% Default settings: RTthresh < 1 min; MZthresh < 0.025 (m/z); FIthresh < 1000
% Later each setting can be set as desired.
[refSet,targetSet,Xr_connIdx,Xt_connIdx, opt]=M2S_matchAll(refFeatures,targetFeatures);
        

% NOTES on setting thresholds:
% Thresholds define the maximum distances that the ref and target features can be in order to be matched.
% You can define "horizontal line" fixed thresholds by using only "intercept" values, and setting "slope" to zero.
% You can define "diagonal line" relative thresholds by using also "slope" values.
% HINT: the function below helps to calculate an intercept and slope in one of the dimensions:
%[RTMZFI_intercept,RTMZFI_slope] = M2S_calculateInterceptSlope() % to define it with the mouse on a figure
%[RTMZFI_intercept,RTMZFI_slope] = M2S_calculateInterceptSlope(startPoint,endPoint) % to insert two points manually

% M2S_matchAll SETTINGS:

% Example settings for SLNEG test data
opt.FIadjustMethod = 'median'; % {'none','median','regression'}
opt.multThresh.RT_intercept = [-0.55,0.15];
opt.multThresh.RT_slope = [0 0];
opt.multThresh.MZ_intercept = [-0.01 0.01]; % m/z units
opt.multThresh.MZ_slope = [-5e-6 5e-6]; % ppm
opt.multThresh.log10FI_intercept = [-1000 1000];
opt.multThresh.log10FI_slope = [0 1];

% Example settings for SLPOS test data
opt.FIadjustMethod = 'median'; % {'median','regression'}
opt.multThresh.RT_intercept = [-0.1,0.1];
opt.multThresh.RT_slope = [0 0];
opt.multThresh.MZ_intercept = [-0.0075 0.015];
opt.multThresh.MZ_slope = [-0.01/2000, -0.01/2000];
opt.multThresh.log10FI_intercept = [-1.4 2.5];
opt.multThresh.log10FI_slope = [0 0];      


% Match all features within defined thresholds

% Define the plot type as:
% No plot: plotType = 0
% Scatter plot: plotType = 1
% Multiple match plot: plotType = 2 (with lines connecting multiple matches containing the same features)

plotType = 1; % No plot = 0;  Scatter plot = 1, or Multiple match plot = 2
[refSet,targetSet,Xr_connIdx,Xt_connIdx,opt]=M2S_matchAll(refFeatures,targetFeatures,opt.multThresh,opt.FIadjustMethod,plotType);

% Obtain properties of current network of matches
[G1,CC] = M2S_infoClusters(refSet,targetSet,1,':');

%%*************************************************************************
%% Procedure part 2: Calculate penalisation scores for each match
%%*************************************************************************

%% Find neighbours, the RT/MZ/log10FI trends, and the residuals - SETTINGS OPTIONAL
% - The objective here is to find the shift trends in RT/MZ/log10FI, and 
% use it to calculate the residuals (distance to those trends).
% - The method to find RT/MZ/FI shift trends uses subsets of features close to  each feature (the "neighbors").
% - The number of neighbors is not critical and can be defined by a percentage of the total of
% features in the refSet (e.g. 0.01) or by a specific number (e.g. 25).
% - The neighbor method to calculate trends can be "cross" (uses RT, MZ, FI)
% or (default) "circle" (uses only RT, MZ). "Circle" does not subtract the trend from FI.
%
% opt.neighbours.nrNeighbors = 0.01; %as a percentage of the number of features in refSet  
% opt.neighbours.nrNeighbors = 21; % as a specific number of neighbours
% opt.calculateResiduals.neighMethod = 'cross'; %{'cross','circle'}

% [Residuals_X,Residuals_trendline] = M2S_calculateResiduals(refSet,targetSet,Xr_connIdx,Xt_connIdx,opt.neighbours.nrNeighbors, opt.calculateResiduals.neighMethod,1)
[Residuals_X,Residuals_trendline] = M2S_calculateResiduals(refSet,targetSet,Xr_connIdx,Xt_connIdx);


%% Adjust the residuals
% - The residuals are in different units and need to be standardised
% to be combined into the penalisation scores for each match.
% (similar to dividing by standard deviation but using double MAD instead)
% by the MAD of the residuals, or by the value at a specific percentile (e.g.75%)
% After this division, the value of the residual at that percentile is 1, so it is easy to visualise the relevance of each of the dimension's residuals (RT/MZ/log10FI).
%
% Uses:
% A: opt.adjustResiduals.residPercentile = NaN % automatic determination using double MAD of Residuals_X. Different value for each dimension. BEST!
% B: opt.adjustResiduals.residPercentile = [0.05,0.004,1.2]; % These selected residuals values become one (residuals from previous figure "Residuals for RT adn MZ")  
% C: opt.adjustResiduals.residPercentile = 80;% Percentile value defined by user, can be between ]0,100] 

%[adjResiduals_X,residPercentile] = M2S_adjustResiduals(refSet,targetSet,Residuals_X,opt.adjustResiduals.residPercentile)
[adjResiduals_X,residPercentile] = M2S_adjustResiduals(refSet,targetSet,Residuals_X);


%% Adjust the weight of each dimension (RT, MZ, log10FI), get penalisation scores
% - It is now possible to set a weight for each of the residuals'
% dimensions. Many times it is desired that log10FI has weight 0.
% - The penalisation scores for each match are the squared root of the sum 
% of squares of the weighted residuals.
%
% opt.weights.W = [1,1,0]; % (Default) equal weight for RT and MZ, not using FI
% opt.weights.W = [20,1,0]; % different weights
% opt.weights.W = [1,1,1]; % equal weight
%
% [penaltyScores] = M2S_defineW_getScores(refSet,targetSet,adjResiduals_X,opt.weights.W,1)
[penaltyScores] = M2S_defineW_getScores(refSet,targetSet,adjResiduals_X);


%% Decide the best of the multiple matches 
% -This recursive method focus in a multiple-match cluster at a time, and 
% finds in each iteration the match with the lowest penalisation score,
% until there are no more selections to do in that cluster.
% - During the calculations you can visualise clusters with multiple matches 
% by giving minNrNodesToPlot a small value, e.g. minNrNodesToPlot = 4.
% 
% plotNetOrNot = 1; % 0 or 1 to plot a global network of all matches colored 
% by the penalisation scores.
%
% [eL,eL_INFO,CC_G1]= M2S_decideBestMatches(refSet, targetSet, Xr_connIdx,Xt_connIdx, penaltyScores, plotNetOrNot, minNrNodesToPlot)
[eL,eL_INFO,CC_G1]= M2S_decideBestMatches(refSet, targetSet, Xr_connIdx,Xt_connIdx, penaltyScores);

%%*************************************************************************
%% Procedure part 3: find false positives (tighten thresholds)
%%*************************************************************************
% - The thresholds set initially are linear, and may not adapt well to
% curved trends. Here we define tighter thresholds that can adapt to curvature.
% - This will delete the matches with largest residuals.
% - Notice the option 'none', which skips this step. Default is 'residuals_mad'
%
% opt.falsePos.methodType = 'scores'; {'none','scores','byBins','trend_mad','residuals_mad'} 
% opt.falsePos.nrMad = 5;
% plotOrNot = 1;
%
% [eL_final, eL_final_INFO] =M2S_findPoorMatches(eL,refSet,targetSet,opt.falsePos.methodType,opt.falsePos.nrMad,plotOrNot);
[eL_final, eL_final_INFO] = M2S_findPoorMatches(eL,refSet,targetSet);

% Summary with number of (multiple matches) discarded, false positive and true positive matches
tableOfMatches = array2table([nansum(isnan(eL_final.notFalsePositives)),nansum(eL_final.notFalsePositives==0),nansum(eL_final.notFalsePositives==1)],'VariableNames',{'DiscardedMatches','PoorMatches','PositiveMatches'});

%%*************************************************************************
%% RESULTS
%%*************************************************************************

% Multiple matches not selected, false positive matches and good matches:
disp(tableOfMatches)

% Table with all results : eL_final
disp(eL_final)
disp(eL_INFO)
 
% The options used in the methods:
disp(opt)

% THE INDICES THAT MATCH refFeatures TO targetFeatures ARE:

refFeatures_idx = eL_final.Xr_connIdx(eL_final.notFalsePositives == 1);
targetFeatures_idx = eL_final.Xt_connIdx(eL_final.notFalsePositives == 1);


%% TO REORDER THE METABOLITE TABLE ACCORDING TO THESE RESULTS (BOTH TO NUMERICAL VALUES AND FEATURE INFO)
% - After we have matched features, we need to select the corresponding
% data from spreadsheets/text files.
% - The following lines produce MZRTstr which are strings with a unique identifier for each metabolomic feature.
% - They also produce vectors of corresponding indices which allow to match the intensity data of reference and target datasets.
% NOTE: features that did not match anything have the index NaN (empty) and can be deleted, as they did not match.

sorted_refIndices_idx = (1:length(refFeatures_idx))';
orderForData_Ref = NaN(size(refFeatures,1),1);
orderForData_Ref(refFeatures_idx) = sorted_refIndices_idx; 
orderForData_Ref_MZRTstr = M2S_createLabelMZRT('REF',refFeatures(:,2),refFeatures(:,1));
matchedRefData_idx = table(orderForData_Ref,orderForData_Ref_MZRTstr,'VariableNames',{'matchedIdxRef','MZRT_ref'});

sorted_targetIndices_idx = (1:length(targetFeatures_idx))';
orderForData_Target = NaN(size(targetFeatures,1),1);
orderForData_Target(targetFeatures_idx) = sorted_targetIndices_idx;
orderForData_Target_MZRTstr = M2S_createLabelMZRT('TARGET',targetFeatures(:,2),targetFeatures(:,1));
matchedTargetData_idx = table(orderForData_Target,orderForData_Target_MZRTstr,'VariableNames',{'matchedIdxRef','MZRT_target'});


%% Write a file with these two vectors:
writetable(matchedRefData_idx,'matchedRefData_idx.xlsx')
writetable(matchedTargetData_idx,'matchedTargetData_idx.xlsx')

% To sort two data files with data (REF and TARGET) according to these indices:
% 0. Run the functions above to save two tables of "matchedData_idx" (ref and target).
% 1. Copy 'matchedRefData_idx' and 'matchedTargetData_idx' to the datasheets (e.g. Excel, csv, etc)
%    of ref and target you want to match, and paste them along the variables
% 2. Sort the variables in ascending order of the NUMERICAL indices (min to max), so the row order becomes 1, 2, 3, etc. The empty indices will be sent to the end. 
% 3. Delete the data for variables with empty indices (these features do not match anything)
% 4. Check there is the same number of variables in both datasets and these
%    have similar (matched) RT and MZ (and FI?). Check that the MZRT_string matches
%    the retention time and m/z of the corresponding features.
% DONE!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SUPPORT FUNCTIONS

[newOpt] = M2S_addOptions(optOld,additionalOpt)
% Objective: add single options to default options without having to define all option values.
% NOTE: all opts are struct variables. 
% Example: additionalOpt.MZ_intercept = 0.05
% [newOpt] = M2S_addOptions(optOld,additionalOpt)

M2S_colorByY_ofSubplot(subplotNr,figH)
% Objective: Color all subplots in figure by y-values of a defined subplot.
% M2S_colorByY_ofSubplot(2,gcf)% colors subplots in current figure
% according to y-value of subplot number 2

[interceptSlope] = M2S_calculateInterceptSlope(startPoint,endPoint)
% Objective: calculate intercept and slope defined by two points.
% Two ways to use:
% 1- By clicking on two points in a plot:
% RTMZFI_slope = M2S_calculateInterceptSlope()
% 2- By inserting the values of the two points ([x1,y1],[x2,y2]) as below:
% RTMZFI_slope = M2S_calculateInterceptSlope([0,0.5],[12,-4.5]);

%% SUPPORT FUNCTIONS FOR OPTIMISATION

[refFeatures_noBigClusters,targetFeatures_noBigClusters] = M2S_deleteLargeClusters(refSet,targetSet,maxFeaturesInCluster,opt);
% Objective: Delete large clusters of multiple matches.
% In some cases the matches contain very large clusters, and it would be
% beneficial/easier to not use those to calculate initial thresholds. This
% function can be used after 'M2S_matchAll', to delete clusters with more
% than e.g. maxFeaturesInCluster = 5. 


[genAlg_Res,optBest] = M2S_genAlg_optimisation(refSet, targetSet, RTdist, MZdist, opt, plotOrNot)
% Objective: optimise initial RT, MZ thresholds.
% Use after running the function M2S_matchAll 
% It runs a genetic algorithm for optimisation.


[G1,CC] = M2S_infoClusters(refSet,targetSet,dimNr,line_style)
% Objective: obtain properties of current network of matches.
% dimNr can be RT=1; MZ=2; FI=3
% example: [G1,CC] = M2S_infoClusters(refSet,targetSet,2,'-')
% disp(CC.freq_clustersWithNedges)



