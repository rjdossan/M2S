%% How to use this function:
%{
% Calculate residuals as usual:
% Match everything:
[refSet,targetSet,Xr_connIdx,Xt_connIdx,opt]=M2S_matchAll(refFeatures,targetFeatures,opt.multThresh,opt.FIadjustMethod,plotType);
[Residuals_X,Residuals_trendline] = M2S_calculateResiduals(refSet,targetSet,Xr_connIdx,Xt_connIdx,opt.neighbours.nrNeighbors, opt.calculateResiduals.neighMethod,opt.pctPointsLoess,plotTypeResiduals)
% To calculate the trendline on RT (dimension 1) on the plot from M2S_matchAll 
[mousePointCoords, predictedPchip] = M2S_manualTrendline(refSet(:,1))
Residuals_trendline(:,1) = predictedPchip;
[Residuals_X,Residuals_trendline] = M2S_interpolateTrendline(refSet,Residuals_trendline,refSet,targetSet)
% Now just adjust the residuals as usual:
opt.adjustResiduals.residPercentile = [0.1,0.0006,25];
[adjResiduals_X,residPercentile] = M2S_adjustResiduals(refSet,targetSet,Residuals_X,opt.adjustResiduals.residPercentile);

%}

function [Residuals_X,Residuals_trendline] = M2S_interpolateTrendline(reduced_refSet,Reduced_trendline,refSet,targetSet)


% loess_IQCimputed = smooth(iQC_runOrder,iQCdata_noLowOutliers_imputed(:,varNr),nrPoints,'loess');% to do loess
% valuesOfBio_inIQC = interp1(iQC_runOrder,loess_IQCimputed,bio_runOrder,'pchip',NaN);
Residuals_X=NaN(size(refSet,1),3);
Residuals_trendline=NaN(size(refSet,1),3);
for v = 1:3
    [~,uniqueIdx]= unique(reduced_refSet(:,v),'stable');
    %Residuals_trendline(:,v) = interp1(reduced_refSet(:,v),Reduced_trendline(:,v),refSet(:,v),'pchip',NaN);
    %Residuals_trendline(:,v) = interp1(reduced_refSet(uniqueIdx,v),Reduced_trendline(uniqueIdx,v),refSet(:,v),'pchip',NaN);
    Residuals_trendline(:,v) = interp1(reduced_refSet(uniqueIdx,v),Reduced_trendline(uniqueIdx,v),refSet(:,v),'pchip','extrap');
    if v<=2
    Residuals_X(:,v) = (targetSet(:,v)-refSet(:,v)) - Residuals_trendline(:,v);
    else
        Residuals_X(:,v) = (log10(targetSet(:,v))-log10(refSet(:,v)))-Residuals_trendline(:,v);
    end
end

%% Sort all dimensions by RT, only to be able to plot lines
M2S_figureH(0.8,0.4)
set(gcf,'name','Final trendlines in each dimension')
for v=1:3
    % Get the reduced sets sorted by reference
    [reduced_refSet_sortedVar, reduced_refSet_sortedVar_idx] = sort(reduced_refSet(:,v));
    Reduced_trendline_sortedVar = Reduced_trendline(reduced_refSet_sortedVar_idx,v);

    % Get the normal sets sorted by reference
    [refSet_sortedVar, refSet_sortedVar_idx] = sort(refSet(:,v));
    Residuals_trendline_sortedVar = Residuals_trendline(refSet_sortedVar_idx,v);

    % Plot
    if v<=2
    subplot(1,3,v), plot(refSet(:,v),targetSet(:,v)-refSet(:,v),'.k'), hold on
    plot(reduced_refSet_sortedVar,Reduced_trendline_sortedVar,'x-b')
    % plot(reduced_refSet(:,v),Reduced_trendline(:,v),'xb')
    hold on, plot(refSet_sortedVar,Residuals_trendline_sortedVar,'-or')
    grid on
    else
        subplot(1,3,v), plot(log10(refSet(:,v)),log10(targetSet(:,v))-log10(refSet(:,v)),'.k'), hold on
        plot(log10(reduced_refSet_sortedVar),(Reduced_trendline_sortedVar),'x-b')
        % plot(reduced_refSet(:,v),Reduced_trendline(:,v),'xb')
        hold on, plot(log10(refSet_sortedVar),(Residuals_trendline_sortedVar),'-or')
    grid on
    end
    axis tight
end
disp('Done!')


M2S_figureH(0.8,0.4)
set(gcf,'name','Final residuals in each dimension')
for v=1:3
    % Plot
    if v<=2
        subplot(1,3,v), plot(refSet(:,v),targetSet(:,v)-refSet(:,v)-Residuals_trendline(:,v),'.k'), hold on
        grid on
    else
        subplot(1,3,v), plot(log10(refSet(:,v)),log10(targetSet(:,v))-log10(refSet(:,v))-Residuals_trendline(:,v),'.k'), hold on       
        grid on
    end
    axis tight
    disp('Done!')
end
disp('Done!')

