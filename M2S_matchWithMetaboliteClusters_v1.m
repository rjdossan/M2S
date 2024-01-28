%% Match FINGER with TILDA

platformType = 'LPOS'
platformType = 'LNEG'


%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load ref FINGER

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\FINGER\IK_analysis\210629MixedModelsIK_BEST\Discussion')
load('WS_220111.mat')
close % 2nd UMAP figure
clearvars -except h_figUMAP plottingData annotationsForPlot groupColorForPlot Xtomatch allPlatforms platformType
% Load updated FINGER annotations
cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\infoOnAllSets\MSannotations')
VarInfo_updated = readtable('230608_FINGER_untargeted_LCMS_VarInfo.xlsx','Sheet','FINGER VarInfo');
% Load the predicted classes using DMA
cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\FINGER\IK_analysis\210629MixedModelsIK_BEST\Discussion\annotation\230612_annotationsForArticle')
VarInfo_predictedClass = tableCellToString(readtable('12-Jun-2023_FINGER_predictedClasses_fromAW1matching.xlsx','sheet','RefMet_ofVarInfo'));
% Load predicted annotations using DMA
cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\FINGER\IK_analysis\210629MixedModelsIK_BEST\Discussion\annotation\230612_annotationsForArticle\230721predictedAnnotationLabels')
tempAnnot = readtable('labelPrediction.xlsx','sheet','VarInfoPredicted');
% Concatenate untargeted predicted annotations
Xtomatch.VarInfo  = tableCellToString([VarInfo_updated,VarInfo_predictedClass(:,2:end),tempAnnot(:,2)]);

% Select only features of the desired platform

idxPlatform = find(contains(Xtomatch.VarInfo.MZRT_str,platformType));
Xtomatch.Data = Xtomatch.Data(:,idxPlatform);
Xtomatch.VarInfo = Xtomatch.VarInfo(idxPlatform,:);


Xref = Xtomatch;
clearvars -except h_figUMAP plottingData Xref platformType

%% Load ref MESA2

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\matchMSfeaturesArticle\matchArticlesFred\testData')

load(['MESA2_S',platformType])
% load MESA2xcms_SLPOS_scaled
Xref.VarInfo.MZRT_str = M2S_createLabelMZRT(['REF_',platformType],Xref.VarInfo.mzmed,Xref.VarInfo.rtmed)
Xref.VarInfo = tableCellToString(Xref.VarInfo);




%% Load target TILDA

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\TILDA\DATA\TILDAxcms_forMatlab')
load(['TILDAxcms_P',platformType,'_scaled_temp'],'X')
X.Data = log(X.Data+1);
X.Data = preprocess_meg(X.Data,4);
Xtarget = X;
Xtarget.VarInfo.MZRT_str = M2S_createLabelMZRT(['T',platformType],Xtarget.VarInfo.mzmed,Xtarget.VarInfo.rtmed);
clearvars -except h_figUMAP plottingData Xref Xtarget platformType

%% Load target Airwave1

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\AIRWAVE\DATA\AIRWAVE1\AIRWAVE1xcms_forMatlab_matched\Xmatched')
load(['Airwave1xcms_S',platformType,'_scaled'],'Xscaled')
% X.Data = log(X.Data+1);
% X.Data = preprocess_meg(X.Data,4);
Xtarget = Xscaled;
Xtarget.VarInfo.MZRT_str = M2S_createLabelMZRT(['T',platformType],Xtarget.VarInfo.mzmed,Xtarget.VarInfo.rtmed);
clearvars -except h_figUMAP plottingData Xref Xtarget platformType

%% Load target MESA phase 2

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\COMBI-BIO\DATA\COMBIBIO2_forMatlab_DC')
load(['MESA2xcms_S',platformType,'_scaled'],'Xscaled')
% X.Data = log(X.Data+1);
% X.Data = preprocess_meg(X.Data,4);
Xtarget = Xscaled;
Xtarget.VarInfo.MZRT_str = M2S_createLabelMZRT(['T',platformType],Xtarget.VarInfo.mzmed,Xtarget.VarInfo.rtmed);
clearvars -except h_figUMAP plottingData Xref Xtarget platformType

%% Load target Rotterdam2

cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\matchMSfeaturesArticle\matchArticlesFred\testData')

load(['ROTTERDAM2_S',platformType])
% load MESA2xcms_SLPOS_scaled
Xtarget.VarInfo.MZRT_str = M2S_createLabelMZRT(['TARGET_',platformType],Xtarget.VarInfo.mzmed,Xtarget.VarInfo.rtmed)
Xtarget.VarInfo = tableCellToString(Xtarget.VarInfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the initial features
refFeatures = [Xref.VarInfo.rtmed,Xref.VarInfo.mzmed,Xref.VarInfo.fimed];
targetFeatures = [Xtarget.VarInfo.rtmed,Xtarget.VarInfo.mzmed,Xtarget.VarInfo.fimed];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot

% Hplot = M2S_plotMZRT_featureSet(featureSet,FIcolorOrNot,marker_size,log10orNot)
figure, subplot(2,1,1), M2S_plotMZRT_featureSet(refFeatures,1,7,1), title('initial ref set')
xlim1 = xlim;
subplot(2,1,2), M2S_plotMZRT_featureSet(targetFeatures,1,7,1), title('initial target set')
xlim2 = xlim;

%% FIND ALL FEATURES OF THE SAME METABOLITE IN EACH DATASET

%% Visualise isotopologue matches to help define thresholds for each dataset 
% NOTE: make tight thresholds

% Test here

refXthresh.RT = 0.01;
refXthresh.MZ = 0.01;
refXthresh.Corr = 0.8;


refXthresh.RT = 0.003;
refXthresh.MZ = 0.004;
refXthresh.Corr = 0.9;


% Define here

refXthresh.RT = 0.001;
refXthresh.MZ = 0.0021;
refXthresh.Corr = 0.85;
[ref_t_rows,ref_r_cols,ref_tr_corr] = DMA_findIsotopologueMatches_v2(Xref,refXthresh,repmat(string(platformType),size(Xref.VarInfo,1),1),1);

targetXthresh.RT = 0.005;
targetXthresh.MZ = 0.0015;
targetXthresh.Corr = 0.95;
[target_t_rows,target_r_cols,target_tr_corr] = DMA_findIsotopologueMatches_v2(Xtarget,targetXthresh,repmat(string(platformType),size(Xtarget.VarInfo,1),1),1);


%% Find all features of the same metabolite

G_metabolite_All_ref = DMA_findAllFeaturesOfMetabolite(Xref, refXthresh, 'undir',1)
G_metabolite_All_target = DMA_findAllFeaturesOfMetabolite(Xtarget, targetXthresh, 'undir',1)

tabulateFormatted(G_metabolite_All_ref.Nodes.metaboliteNr)
figure, histogram(G_metabolite_All_ref.Nodes.metaboliteNr,'BinMethod','integers')
tabulateFormatted(G_metabolite_All_target.Nodes.metaboliteNr)

% Find idx of network nodes in original set
G_metabolite_All_ref.Nodes.idxInVarInfo = M2S_find_idxInReference(G_metabolite_All_ref.Nodes.Name,Xref.VarInfo.MZRT_str)
G_metabolite_All_target.Nodes.idxInVarInfo = M2S_find_idxInReference(G_metabolite_All_target.Nodes.Name,Xtarget.VarInfo.MZRT_str)

%% match all features in a metabolite

if datasetName == "LNEG"
    opt.FIadjustMethod = 'median'; % {'median','regression'}
    opt.multThresh.RT_intercept = [0,0.15];
    opt.multThresh.RT_slope = [0 0];
    opt.multThresh.MZ_intercept = [-0.0075 0.003];
    opt.multThresh.MZ_slope = [0 0];
    opt.multThresh.log10FI_intercept = [-2 2];
    opt.multThresh.log10FI_slope = [0 0]; 
    
elseif datasetName == "LPOS"
    opt.FIadjustMethod = 'median'; % {'median','regression'}
    opt.multThresh.RT_intercept = [-0.2,0];
    opt.multThresh.RT_slope = [0 0];
    opt.multThresh.MZ_intercept = [-0.009 0.006];
    opt.multThresh.MZ_slope = [0 0];
    opt.multThresh.log10FI_intercept = [-20 20];
    opt.multThresh.log10FI_slope = [0 0];    
end

% MESA2 vs ROTTERDAM LPOS
opt.FIadjustMethod = 'median'; % {'median','regression'}
opt.multThresh.RT_intercept = [-0.1,0.1];
opt.multThresh.RT_slope = [0 0];
opt.multThresh.MZ_intercept = [-0.0075 0.015];
opt.multThresh.MZ_slope = [-0.01/2000, -0.01/2000];
opt.multThresh.log10FI_intercept = [-1.4 0.75];
opt.multThresh.log10FI_slope = [0 0];    

% MESA2 vs ROTTERDAM LNEG
opt.FIadjustMethod = 'median'; % {'median','regression'}
opt.multThresh.RT_intercept = [-1,1];
opt.multThresh.RT_slope = [0 0];
opt.multThresh.MZ_intercept = [-0.01 0.01];
opt.multThresh.MZ_slope = [0 0];
opt.multThresh.log10FI_intercept = [-2 2];
opt.multThresh.log10FI_slope = [0 0]; 

% Match features within defined thresholds (all features, published method)
plotType = 3; 
[refSet,targetSet,Xr_connIdx,Xt_connIdx,opt]=M2S_matchAll(refFeatures,targetFeatures,opt.multThresh,opt.FIadjustMethod,plotType);


% Match features within defined thresholds (only the ones in metabolites)
% [refSet_ap,targetSet_ap,Xr_connIdx_local_ap,Xt_connIdx_local_ap,opt_local_ap]=M2S_matchAll(refFeatures(G_metabolite_All_ref.Nodes.idxInVarInfo,:),targetFeatures(G_metabolite_All_target.Nodes.idxInVarInfo,:));
plotType = 3; 
[refSet_withMultiples,targetSet_withMultiples,Xr_connIdx_withMultiples_local,Xt_connIdx_withMultiples_local,opt]=M2S_matchAll(refFeatures(G_metabolite_All_ref.Nodes.idxInVarInfo,:),targetFeatures(G_metabolite_All_target.Nodes.idxInVarInfo,:),opt.multThresh,opt.FIadjustMethod,plotType);
figure, M2S_plotDelta_matchedSets(refSet_withMultiples,targetSet_withMultiples,'.k',2)

% Evalutate with the correlated features
[ResCFSM] = match_using_correlatedFeatures(X1,X2,Xr_connIdx,Xt_connIdx,scoreToColour)



%% 
%% Find another way to select only relevant features after multiple matching matches
figure, M2S_plotDelta_matchedSets(refSet,targetSet,'.k',2)
hold on, M2S_plotDelta_matchedSets(refSet_withMultiples,targetSet_withMultiples,'or',2)
% Select only large peaks
medianUniqueFIref = median(unique(refSet(:,3)))
refSet_largePeaks_localIdx = find(refSet(:,3)>medianUniqueFIref)

medianUniqueFItarget = median(unique(targetSet(:,3)))
targetSet_largePeaks_localIdx = find(targetSet(:,3)>medianUniqueFItarget)

refTargetSet_largePeaks_localIdx = intersect(refSet_largePeaks_localIdx, targetSet_largePeaks_localIdx)

refSet_largePeaks_globalIdx = Xr_connIdx(refTargetSet_largePeaks_localIdx)
targetSet_largePeaks_globalIdx = Xt_connIdx(refTargetSet_largePeaks_localIdx)

hold on, M2S_plotDelta_matchedSets(refSet(refTargetSet_largePeaks_localIdx,:),targetSet(refTargetSet_largePeaks_localIdx,:),'or',2)


% Select only unique matches

apR = tabulateFormatted(string(Xr_connIdx))
refSet_uniqueMatches_globalIdx = str2double(table2array(apR(table2array(apR(:,2) == 1),1)))
refSet_uniqueMatches_localIdx = find(M2S_find_idxInReference(Xr_connIdx, refSet_uniqueMatches_globalIdx))

apT = tabulateFormatted(string(Xt_connIdx))
targetSet_uniqueMatches_globalIdx = str2double(table2array(apT(table2array(apT(:,2) == 1),1)))
targetSet_uniqueMatches_localIdx = find(M2S_find_idxInReference(Xt_connIdx, targetSet_uniqueMatches_globalIdx))

refTargetSet_uniqueMatches_localIdx = intersect(refSet_uniqueMatches_localIdx, targetSet_uniqueMatches_localIdx)






% Find the global idx and MZRTstr
Xr_connIdx_withMultiples = G_metabolite_All_ref.Nodes.idxInVarInfo(Xr_connIdx_withMultiples_local);
Xt_connIdx_withMultiples = G_metabolite_All_target.Nodes.idxInVarInfo(Xt_connIdx_withMultiples_local);
MZRTstr_withMultiplesRef = Xref.VarInfo.MZRT_str(Xr_connIdx_withMultiples);
MZRTstr_withMultiplesTarget = Xtarget.VarInfo.MZRT_str(Xt_connIdx_withMultiples);

% Find the metabolite number (ConnectedComponent) of each matched feature
CC_connRef = G_metabolite_All_ref.Nodes.metaboliteNr(M2S_find_idxInReference(MZRTstr_withMultiplesRef, string(G_metabolite_All_ref.Nodes.Name)))
CC_connTarget = G_metabolite_All_target.Nodes.metaboliteNr(M2S_find_idxInReference(MZRTstr_withMultiplesTarget, string(G_metabolite_All_target.Nodes.Name)))

matchedMZRT_str_Ref = [];
matchedMZRT_str_Target = [];
unique_CC_connRef = unique(CC_connRef)
for ccidx = 1:length(unique_CC_connRef)
    tempIdx = find(CC_connRef == unique_CC_connRef(ccidx));
    Tab_freqCCtarget = tabulateFormatted(CC_connTarget(tempIdx));
    if Tab_freqCCtarget.Count(1) > 1
        selectedTargetCluster = str2double(Tab_freqCCtarget.Value(1));
        matches_idxInVarInfo = find(CC_connRef == unique_CC_connRef(ccidx) & CC_connTarget == selectedTargetCluster)
        [MZRTstr_withMultiplesRef(matches_idxInVarInfo), MZRTstr_withMultiplesTarget(matches_idxInVarInfo)];
        matchedMZRT_str_Ref = [matchedMZRT_str_Ref;MZRTstr_withMultiplesRef(matches_idxInVarInfo)];
        matchedMZRT_str_Target = [matchedMZRT_str_Target; MZRTstr_withMultiplesTarget(matches_idxInVarInfo)];
    end
end
Xr_connIdx_withMultiplesSelected = M2S_find_idxInReference(matchedMZRT_str_Ref,Xref.VarInfo.MZRT_str)
Xt_connIdx_withMultiplesSelected = M2S_find_idxInReference(matchedMZRT_str_Target,Xtarget.VarInfo.MZRT_str)

% Plot results
% [refSet,targetSet,Xr_connIdx,Xt_connIdx, opt]=M2S_matchAll(refFeatures,targetFeatures);
     
refSet_withMultiplesSelected_rawFI = refFeatures(Xr_connIdx_withMultiplesSelected,:);
targetSet_withMultiplesSelected_rawFI = targetFeatures(Xt_connIdx_withMultiplesSelected,:);

refSet_rawFI = refFeatures(Xr_connIdx,:);
targetSet_rawFI = targetFeatures(Xt_connIdx,:);

% Plot

figure
plotRow = 2;
Xsymbol = '.k';  % y-axis is e.g. RTdist (1) or both RTdist and RT_target (2)
M2S_plotDelta_matchedSets(refSet_rawFI,targetSet_rawFI,Xsymbol,plotRow)
Xsymbol = 'or';  % y-axis is e.g. RTdist (1) or both RTdist and RT_target (2)
M2S_plotDelta_matchedSets(refSet_withMultiplesSelected_rawFI,targetSet_withMultiplesSelected_rawFI,Xsymbol,plotRow)


% Rematch to harmonise the FI
plotType = 3; 
[refSet_withMultiplesSelected,targetSet_withMultiplesSelected,~,~,~]=M2S_matchAll(refFeatures(Xr_connIdx_withMultiplesSelected,:),targetFeatures(Xt_connIdx_withMultiplesSelected,:),opt.multThresh,opt.FIadjustMethod,plotType);

%% Create penalisation scores using function to interpolate (rloess + pchip)

% 1. This first one is only on features from "metabolites", to calculate a threshold
opt.neighbours.nrNeighbors = 0.2;
[Residuals_X_forThreshold,~] = M2S_interpolate(refSet_withMultiplesSelected,targetSet_withMultiplesSelected,refSet_withMultiplesSelected,targetSet_withMultiplesSelected, opt.neighbours.nrNeighbors)

residualsMAD_forThreshold = mad(Residuals_X_forThreshold)
% residualsSTD = std(Residuals_X_forThreshold)

% 2. This is the one to obtain the residuals and inter-dataset shift trend lines

% NOTE: First change the values of refSet outside of range of refSet_withMultiplesSelected
% to the minimum or maximum value of refSet_withMultiplesSelected. This is
% because otherwise the values outside the range can give weird residuals.
 
refSet_inRange = refSet;
for d=1:3
    refSet_inRange(refSet_inRange(:,d)<min(refSet_withMultiplesSelected(:,d)),d) = min(refSet_withMultiplesSelected(:,d));
    refSet_inRange(refSet_inRange(:,d)>max(refSet_withMultiplesSelected(:,d)),d) = max(refSet_withMultiplesSelected(:,d));
end

opt.neighbours.nrNeighbors = 0.2;
[Residuals_X,Residuals_trendline] = M2S_interpolate(refSet_withMultiplesSelected,targetSet_withMultiplesSelected,refSet_inRange, targetSet, opt.neighbours.nrNeighbors)
residualsMAD = mad(Residuals_X)

% Calculate the ratio between MAD of selected or all features, and its mean value (only for RT and MZ)
ResidualsMAD_ratio = residualsMAD_forThreshold./residualsMAD;
dimensionsToAverageMAD = 1:2; % without FI = 1:2, with FI = 1:3
ResidualsMAD_ratio_mean = sqrt(mean((ResidualsMAD_ratio(dimensionsToAverageMAD)).^2));

% Normalise the residuals at 1*MAD (calculated only from features in metabolites)
[adjResiduals_X,residPercentile] = M2S_adjustResiduals(refSet, targetSet,Residuals_X, residualsMAD_forThreshold);

% Create penalty scores
opt.weights.W = [1,1,0]; % equal weight to RT, MZ, no weight to FI
[penaltyScores] = M2S_defineW_getScores(refSet, targetSet,adjResiduals_X,opt.weights.W,1); 

% Decide best matches
[eL,eL_INFO,CC_G1]= M2S_decideBestMatches(refSet, targetSet, Xr_connIdx,Xt_connIdx, penaltyScores);

% Find poor matches using the scores method
opt.falsePos.methodType = 'scores'; %{'none','scores','byBins','trend_mad','residuals_mad'} 
opt.falsePos.nrMad = 3;
plotOrNot = 1;
[eL_final, eL_final_INFO] = M2S_findPoorMatches(eL,refSet,targetSet,opt.falsePos.methodType, ResidualsMAD_ratio_mean * opt.falsePos.nrMad,plotOrNot);

tableOfMatches = array2table([nansum(isnan(eL_final.notFalsePositives)),nansum(eL_final.notFalsePositives==0),nansum(eL_final.notFalsePositives==1)],'VariableNames',{'DiscardedMatches','PoorMatches','PositiveMatches'})


% THE INDICES THAT MATCH refFeatures TO targetFeatures ARE:

refFeatures_idx = eL_final.Xr_connIdx(eL_final.notFalsePositives == 1);
targetFeatures_idx = eL_final.Xt_connIdx(eL_final.notFalsePositives == 1);

% Prepare the dataset

tempTargetVarInfo = createEmptyTable(Xref.VarInfo,size(Xtarget.VarInfo,1));
tempTargetVarInfo(targetFeatures_idx,:) = Xref.VarInfo(refFeatures_idx,:);
tempTargetVarInfo.Properties.VariableNames = strcat("FINGER_",string(tempTargetVarInfo.Properties.VariableNames))
TILDAvarinfo_annotatedByFINGER = tableCellToString([Xtarget.VarInfo,tempTargetVarInfo]);
% TILDA_annotatedByFINGER = Xtarget.VarInfo;


% Count annotations
% From FINGER/AIRWAVE 
sum(TILDAvarinfo_annotatedByFINGER.FINGER_AbbreviatedAnnotation~="")
% Predicted by DMA
sum(TILDAvarinfo_annotatedByFINGER.FINGER_AbbreviatedAnnotationFinal~="")


cd('C:\Users\rjdossan\OneDrive - Imperial College London\Imperial projects\TILDA\DATA\TILDAxcms_forMatlab')
writetable(TILDAvarinfo_annotatedByFINGER, 'TILDAvarinfo_annotatedByFINGER.xlsx','Sheet',platformType)




