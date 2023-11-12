%% This is a function to reduce the number of peaks to the most intense.
% It divides refFeatures and targetFeatures per bins in RT and MZ, then 
% keeps the most intense peak in each bin.
% nBinsRef will find that number of bins both in RT and in MZ separately,
% then will consider the most intense peak in each of these bins and merge
% the ones from RT and mz. The final number of bins/features obtained in the 
% end is in general larger than the number defined.
% [refFeatures_reduced,targetFeatures_reduced, refFeaturesBin_maxFI_idx, targetFeaturesBin_maxFI_idx] = M2S_reduceSet(refFeatures, targetFeatures, nBinsRef, nBinsTarget)

function [refFeatures_reduced,targetFeatures_reduced, refFeaturesBin_maxFI_idx, targetFeaturesBin_maxFI_idx] = M2S_reduceSet(refFeatures, targetFeatures, nBinsRef, nBinsTarget)


% nBinsRef = 1000 % n bins refFeatures (always > 1)
% nBinsTarget = 500 % n bins in targetFeatures

% Find the bin for each match set



%% Select the feature with max FI in each bin in Reference
refFeaturesBin_maxFI_idx = [];
refIdxVec = (1:size(refFeatures,1))';
c = 1;
% For Retention time
binPer_refFeature = discretize(refFeatures(:,1),nBinsRef);
for binNr = 1:max(binPer_refFeature)
    binIdx = find(binPer_refFeature == binNr);
    if ~isempty(binIdx)
         [~,localMaxIdx] = max(refFeatures(binIdx,3));
         refFeaturesBin_maxFI_idx(c) = refIdxVec(binIdx(localMaxIdx));
         c = c+1;
    end
end

% For m/z
binPer_refFeature = discretize(refFeatures(:,2),nBinsRef);
for binNr = 1:max(binPer_refFeature)
    binIdx = find(binPer_refFeature == binNr);
    if ~isempty(binIdx)
         [~,localMaxIdx] = max(refFeatures(binIdx,3));
         refFeaturesBin_maxFI_idx(c) = refIdxVec(binIdx(localMaxIdx));
         c = c+1;
    end
end
% Collect indices of all reference features chosen
refFeaturesBin_maxFI_idx = unique(refFeaturesBin_maxFI_idx','stable');

%% Select the feature with max FI in each bin in Target
targetFeaturesBin_maxFI_idx = [];
targetIdxVec = (1:size(targetFeatures,1))';
c = 1;
% For Retention time
binPer_targetFeature = discretize(targetFeatures(:,1),nBinsTarget);
for binNr = 1:max(binPer_targetFeature)
    binIdx = find(binPer_targetFeature == binNr);
    if ~isempty(binIdx)
         [~,localMaxIdx] = max(targetFeatures(binIdx,3));
         targetFeaturesBin_maxFI_idx(c) = targetIdxVec(binIdx(localMaxIdx));
         c = c+1;
    end
end

% For m/z
binPer_targetFeature = discretize(targetFeatures(:,2),nBinsTarget);
for binNr = 1:max(binPer_targetFeature)
    binIdx = find(binPer_targetFeature == binNr);
    if ~isempty(binIdx)
         [~,localMaxIdx] = max(targetFeatures(binIdx,3));
         targetFeaturesBin_maxFI_idx(c) = targetIdxVec(binIdx(localMaxIdx));
         c = c+1;
    end
end

% Collect indices of all target features chosen
targetFeaturesBin_maxFI_idx = unique(targetFeaturesBin_maxFI_idx','stable');

% Collect the reduced Feature sets 
refFeatures_reduced = refFeatures(refFeaturesBin_maxFI_idx,:);
targetFeatures_reduced = targetFeatures(targetFeaturesBin_maxFI_idx,:);



