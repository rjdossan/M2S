% M2S_threshSelection
% This function is used in the context of matching features between untargeted metabolomics datasets.
% It uses a concept similar to the double MAD.
%

function [S_result,percentileThresh] = M2S_threshSelection(S, methodType, pAlpha_or_nrMad, plotOrNot)

% NOTE: if method is "mad", pAlpha_or_nrMad is nrMad (e.g. 3).
% S = RTMZdist_medNeigh_to_TRpair_scaled(eL.is_Best==1);

plotOrNotIterations = 0; % This is in case one wants to plot all the iterations to visualise Method 1.

if nargin ==1
    methodType = 'mad';
    pAlpha_or_nrMad = 3;
    plotOrNot = 1;

elseif nargin == 2
%     if strcmp(methodType,'outliers')
%         pAlpha_or_nrMad = 0.01;
    if strcmp(methodType,'mad')
        pAlpha_or_nrMad = 3;
    else
        disp('The selected methodType does not exist')
    end
    plotOrNot = 1;

elseif  nargin == 3
        plotOrNot = 1;
end

%% ******************************************************************
if strcmp(methodType,'mad')
% Leys, C., et al., Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median, Journal of Experimental Social Psychology, Volume 49, Issue 4, July 2013, pp. 764-766. *
% Rousseeuw, P.J. and Croux C. (1993) Alternatives to the Median Absolute Deviation, Journal of the American Statistical Association, December 1993, pp. 1273-1283.    
    
    % determine the double MAD
    percentForNbins=0.1;
    % pAlpha_or_nrMad = 3;
    % high_S_idx = (1:length(S))'; %find(S>nanmedian(S));
    % high_S = S(high_S_idx);% 
    y_S = S + 0.05 *range(S) .* rand(size(S));% only for the plot

    
    % OUTPUT
    S_limit = (median(S) + mad(S)*pAlpha_or_nrMad);
    S_used_insideOfLimit_idx = find(S< S_limit);
    S_allOutOfLimit_idx = find(S> S_limit);
    S_used_outOfLimit_idx = NaN;
    
    if plotOrNot == 1
        
        figure('Position',[1   445   640   500]) , movegui(gcf,'center')
        subplot(2,1,1)
        plot(S,y_S,'.k'), hold on
        M2S_plotaxes('k',[median(S),NaN]);
        M2S_plotaxes('r',[S_limit,NaN]);
        axis tight, grid on
        [hhh_Values, hhh_BinEdges] = histcounts(S,round(percentForNbins*length(S)),'Normalization','pdf'); 
        hhh_Values(isinf(hhh_Values)) = NaN; hhh_Values = (max(ylim)/nanmax(hhh_Values)) * hhh_Values;
        hhh_BinCentres = hhh_BinEdges(2:end) - 0.5*(hhh_BinEdges(2)-hhh_BinEdges(1));
        bar(hhh_BinCentres,hhh_Values,'b');
        title('MAD method: Scores distribution, median (black line) and threshold (red line)')
        xlabel('Scores'),ylabel('Frequency')
        subplot(2,1,2)
        %figure('Position',[1    97   640   262]), 
        plot(S_used_insideOfLimit_idx,S(S_used_insideOfLimit_idx),'.k'), hold on
        plot(S_allOutOfLimit_idx,S(S_allOutOfLimit_idx),'.r'), grid on, axis tight
        plot((xlim)',S_limit*ones(2,1),'-r')
        title('MAD method: Scores inside limit in black, outside in red')
        xlabel('Variable index'),ylabel('Scores')
    end

end

S_result.S_used_insideOfLimit_idx = S_used_insideOfLimit_idx;
S_result.S_allOutOfLimit_idx = S_allOutOfLimit_idx;
S_result.S_used_outOfLimit_idx = S_used_outOfLimit_idx;
S_result.S_limit = S_limit;% THE THRESHOLD IN SCORES UNITS
percentileThresh = length(S_result.S_used_insideOfLimit_idx)/(length(S_result.S_used_insideOfLimit_idx)+length(S_result.S_allOutOfLimit_idx));
S_result.options.methodType = methodType;
S_result.options.pAlpha_or_nrMad = pAlpha_or_nrMad;
