% M2S_threshSelection
% This function is used in the context of matching features between untargeted metabolomics datasets.
% It is used to automatically select a threshold for penalty scores
% (all positive) by comparison to a chi square distribution. Notice that
% it assumes that the degree of freedom of the chi distribution is 1.
% The difficulty in comparing the distribution of scores with the expected
% values of chi square is that there are outliers. Thus this function
% creates multiple subsets from the scores, from a constant minimum value to a
% maximum value that varies. It then compares each of those datasets with
% the expected values of chi square and selects the closest one.
% NOTE: the e.g. chiCriticalVal == 15 is arbitrarily selected to represent the last possible
% value of the chi distribution (minimum probability, close to zero). Its chiSquare cdf is p = 0.0001075;
%
% Previously it was   function [S_result] = autothreshSelectionautoThreshSelectionChi_v2(S, methodType, pAlpha_or_nrMad, plotOrNot)

function [S_result,percentileThresh] = M2S_threshSelection(S, methodType, pAlpha_or_nrMad, plotOrNot)

% NOTE: if method is "outliers", pAlpha_or_nrMad is pAlpha (e.g. 0.01). 
%       if method is "mad", pAlpha_or_nrMad is nrMad (e.g. 3).
% S = RTMZdist_medNeigh_to_TRpair_scaled(eL.is_Best==1);
% S =  chi2rnd(2,4923,1); % to test with chi square
% distribution and df==1
% pAlpha = 0.01;
plotOrNotIterations = 0; % This is in case one wants to plot all the iterations to visualise Method 1.

if nargin ==1
    methodType = 'mad';
    pAlpha_or_nrMad = 3;
    plotOrNot = 1;

elseif nargin == 2
    if strcmp(methodType,'outliers')
        pAlpha_or_nrMad = 0.01;
    elseif strcmp(methodType,'mad')
        pAlpha_or_nrMad = 3;
    else
        disp('The selected methodType does not exist')
    end
    plotOrNot = 1;

elseif  nargin == 3
        plotOrNot = 1;
end


%% METHOD 1
% [h,p] = chi2gof(r,'cdf',@(x)chi2cdf(x,2),'nparams',1);

if strcmp(methodType,'outliers') % NOT USED!
    % Initialisations
    n_testPoints = 20;
    df = 1; % the degree of freedom of the chi square distribution is set to 1
    percentForNbins = 0.1;
    maxChiCriticalValue = 15;

    S_areaLimits = (linspace(0,max(S),n_testPoints+1))';
    S_areaLimits = S_areaLimits(2:end);% do not choose 0

    chiT=struct; % chiTheoretical values
    h=struct; % structure to collect the scores results
    relSSQ = NaN(length(S_areaLimits),1);% the relative error

    if plotOrNotIterations == 1
        figure('Position',[94          74        1700         888])
    end
    % for each test set
    for f = 1:length(S_areaLimits)
        %% first select a set of S__sorted
        h(f).S_used_initial_idx = find(S<=S_areaLimits(f));
        h(f).S_notUsed_initial_idx = find(S>S_areaLimits(f));

        h(f).S_used = S(h(f).S_used_initial_idx);
        h(f).NumBins = round(percentForNbins*length(h(f).S_used));

        %% then get the theoretical pdf for the critical values of the chisquare distribution
        chiT(f).chiCriticalVals = (linspace(0,maxChiCriticalValue,h(f).NumBins))';
        chiT(f).chi2valsPdf = pdf('Chisquare',chiT(f).chiCriticalVals,df);
        % Normalize: get values of scores in same magnitude as chi2 critical values
        h(f).S_used_norm = (max(chiT(f).chiCriticalVals)/max(h(f).S_used)) * h(f).S_used; 

        [temp_Values,temp_BinEdges] = histcounts(h(f).S_used_norm,round(percentForNbins*length(h(f).S_used_norm)),'Normalization','pdf');
        h(f).BinUpperLimits = temp_BinEdges(2:end)';
        h(f).Values = temp_Values';
        if df==1 % in this case the theoretical chi2 first value is infinite
            chiT(f).chi2valsPdf(1) = h(f).Values(1);
        end
        h(f).normFactor = max(chiT(f).chiCriticalVals)/max(h(f).S_used);
        % Calculate the error for this test set
        relSSQ(f,1) = sum((h(f).Values - chiT(f).chi2valsPdf).^2) ./ sum((chiT(f).chi2valsPdf).^2);
        if plotOrNotIterations == 1
            plot(h(f).BinUpperLimits,h(f).Values,'.k','MarkerSize',8), axis tight
            S_notUsed_norm = (h(f).normFactor) * S(h(f).S_notUsed_initial_idx); % get values of scores in same magnitude as chi2 critical values
            hold on, plot(S_notUsed_norm,zeros(size(S_notUsed_norm)),'xb')% plot the values that were not selected in this iteration
            plot(chiT(f).chiCriticalVals,chiT(f).chi2valsPdf,'.-r'), axis tight; title(['test set ',num2str(f)]); grid on, hold off
            pause; clf
        end
    end    


    % CONCLUSIONS FOR METHOD 1
    [~,bestTestSet] = min(relSSQ);

    % Additional results
    diffTo_pAlpha = cdf('Chisquare', linspace(0,maxChiCriticalValue,h(bestTestSet).NumBins), df) - (1-pAlpha_or_nrMad);
    [diffVal, diffValMin_idx] = min(abs(diffTo_pAlpha));

    % OUTPUT
    S_limit = (h(bestTestSet).BinUpperLimits(diffValMin_idx)/h(bestTestSet).normFactor);% the Scores threshold limit
    S_used_insideOfLimit_idx = find(S <= S_limit);% good ones
    S_allOutOfLimit_idx = find(S > S_limit); % all out of limit
    S_used_outOfLimit_idx = setdiff(S_allOutOfLimit_idx,h(bestTestSet).S_notUsed_initial_idx); % used out of limit

    % PLOTS FOR METHOD 1
    if plotOrNot == 1
        % Plot the error for all the iterations. The minimum the better
        %figure('Position',[1   771   640   224]), 
        figure('Position',[1   1  739   810]), movegui(gcf,'center')
        subplot(3,1,1), plot(relSSQ,'o-k')
        hold on, plot(bestTestSet,relSSQ(bestTestSet),'or'), axis tight, grid on, 
        title('Outliers method: Relative SS error for each iteration')
        xlabel('Iteration number'),ylabel('Relative SS error')

        % Plot the theoretical chi2 and the observed values adjusted as chi2
%         figure('Position',[2   445   640   241]), 
        subplot(3,1,2), title('Outliers method: Best adjusted observed values vs theoretical chi2 values')
        plot(h(bestTestSet).BinUpperLimits,h(bestTestSet).Values,'.k','MarkerSize',8), axis tight, grid on
        hold on, plot(chiT(bestTestSet).chiCriticalVals,chiT(bestTestSet).chi2valsPdf,'.-r')
        M2S_plotaxes('g',[chiT(bestTestSet).chiCriticalVals(diffValMin_idx),NaN])% critical value obtained from theoretical chi2
        %M2S_plotaxes('r',[h(bestTestSet).BinUpperLimits(diffValMin_idx),NaN]); % critical value obtained from histogram
        M2S_plotaxes('k',[(df*(1 - (2/(9*df)))^3),0]);% chi2 median
        xlabel('Theoretical chi2 critical'),ylabel('pdf of observed values')
        
        % get the local idx of scores outside the alpha p value limit
        %figure('Position',[1    97   640   262])
        subplot(3,1,3), 
        plot(h(bestTestSet).S_used_initial_idx,S(h(bestTestSet).S_used_initial_idx),'.k'); hold on % used
        plot(h(bestTestSet).S_notUsed_initial_idx,S(h(bestTestSet).S_notUsed_initial_idx),'.b') % not used
        plot(S_used_outOfLimit_idx,S(S_used_outOfLimit_idx),'or')% used but outside of limit
        grid on, axis tight; 
        plot((xlim)',S_limit*ones(2,1),'-r')
        title('Outliers method: Scores inside limit in black, outside in red, chi2 outliers in blue')       
        xlabel('Variable index'),ylabel('Scores')
    end

%% ******************************************************************
elseif strcmp(methodType,'mad')
% Leys, C., et al., Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median, Journal of Experimental Social Psychology, Volume 49, Issue 4, July 2013, pp. 764-766. *
% Rousseeuw, P.J. and Croux C. (1993) Alternatives to the Median Absolute Deviation, Journal of the American Statistical Association, December 1993, pp. 1273-1283.    
    
%% METHOD 2
    % determine the double MAD
    percentForNbins=0.1;
    %pAlpha_or_nrMad = 3;
    high_S_idx = find(S>nanmedian(S));
    high_S = S(high_S_idx);
    y_S = S + 0.05 *range(S) .* rand(size(S));% only for the plot

    
    % OUTPUT
    S_limit = (median(high_S) + mad(high_S)*pAlpha_or_nrMad);
    S_used_insideOfLimit_idx = find(S< S_limit);
    S_allOutOfLimit_idx = find(S> S_limit);
    S_used_outOfLimit_idx = NaN;
    
    if plotOrNot == 1
        
        figure('Position',[1   445   640   500]) , movegui(gcf,'center')
        subplot(2,1,1)
        plot(S,y_S,'.k'), hold on
        M2S_plotaxes('k',[median(high_S),NaN]);
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
