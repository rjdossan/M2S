%%  [Residuals_X,Residuals_trendline] = M2S_interpolate(refSet,targetSet,refSetToPredict, targetSetToPredict, nPoints)
% This function calculates loess on the first unique points of x, 
% % then predicts all x points using pchip (polynomial cubic interpolation)
% refSet and targetSet may be 3 columns (RT, MZ, FI not log10) or a single
% column with the desired vectors to interpolate.
% nPoints: optional (default is 0.2) (for rloess) can be between 0-1 (percent) or a number of points >1, e.g. 20 points

function [Residuals_X,Residuals_trendline] = M2S_interpolate(refSet,targetSet,refSetToPredict, targetSetToPredict, nPoints)

if nargin == 4
    nPoints = 0.2;
end

if size(refSet,2)>1 % calculate in the three dimensions, RT, MZ, FI
    refSet(:,3) = log10(refSet(:,3));
    targetSet(:,3) = log10(targetSet(:,3));
    refSetToPredict(:,3) = log10(refSetToPredict(:,3));
    targetSetToPredict(:,3) = log10(targetSetToPredict(:,3));

    for d=1:3 
        [~,currentDim_uniqueIdx] = unique(refSet(:,d));
        current_loess = smooth(refSet(currentDim_uniqueIdx,d), targetSet(currentDim_uniqueIdx,d),nPoints,'rloess');% to do loess
        predictedPchip(:,d) = pchip(refSet(currentDim_uniqueIdx,d),current_loess ,refSetToPredict(:,d));
    end
        
else % calculate only on one dimension, actually on the values entered
    [~,currentDim_uniqueIdx] = unique(refSet);
    current_loess = smooth(refSet(currentDim_uniqueIdx), targetSet(currentDim_uniqueIdx),nPoints,'rloess');% to do loess
    predictedPchip = pchip(refSet(currentDim_uniqueIdx),current_loess ,refSetToPredict);
end

Residuals_trendline = predictedPchip;
Residuals_X = targetSetToPredict - Residuals_trendline;

%% PLOT
 %% Figure with the trends for RT, MZ, FI
if size(refSet,2)>1
    M2S_figureH(0.65,0.35); set(gcf,'Name','Trends for RT, MZ and log10FI')
    subplot(1,3,1),
    plot(refSetToPredict(:,1),targetSetToPredict(:,1) - refSetToPredict(:,1),'.k'), axis tight, grid on
    y_lim=ylim; x_lim = xlim; hold on; plot(refSetToPredict(:,1),Residuals_trendline(:,1)-refSetToPredict(:,1) ,'or','MarkerSize',4), grid on
    ylim(y_lim); xlim(x_lim); xlabel('RT reference (Minutes)'), ylabel('RTdist (Minutes)')
    subplot(1,3,2)
    plot(refSetToPredict(:,2),targetSetToPredict(:,2) - refSetToPredict(:,2),'.k'), axis tight, grid on
    y_lim=ylim; x_lim = xlim; hold on; plot(refSetToPredict(:,2),Residuals_trendline(:,2)-refSetToPredict(:,2),'or','MarkerSize',4), grid on
    ylim(y_lim); xlim(x_lim); xlabel('MZ reference (m/z units)'), ylabel('MZdist (m/z units)')
    subplot(1,3,3)
    plot(refSetToPredict(:,3),targetSetToPredict(:,3) - refSetToPredict(:,3),'.k'), axis tight, grid on
    y_lim=ylim; x_lim = xlim; hold on; plot(refSetToPredict(:,3),Residuals_trendline(:,3)-refSetToPredict(:,3),'or','MarkerSize',4), grid on
    ylim(y_lim); xlim(x_lim); xlabel('log10FI reference'), ylabel('log10FIdist')
    drawnow  
    
   
    %% Figure with the Residuals_X
    M2S_figureH(0.45,0.35); set(gcf,'Name','Residuals for RT, MZ and log10FI')
    subplot(1,3,1),
    plot(refSetToPredict(:,1),Residuals_X(:,1),'.k'), axis tight, hold on, xlim1 = xlim; grid on
    %plot(xlim',[0;0],'-k')
    xlabel('RT of reference feature'), ylabel('RT residuals') 
    subplot(1,3,2),plot(refSetToPredict(:,2),Residuals_X(:,2),'.k'), axis tight, hold on, xlim1 = xlim; grid on
    %plot(xlim',[0;0],'-k')
    xlabel('MZ of reference feature'), ylabel('MZ residuals')
    subplot(1,3,3),plot(refSetToPredict(:,3),Residuals_X(:,3),'.k'), axis tight, hold on, xlim1 = xlim; grid on
    %plot(xlim',[0;0],'-k')
    xlabel('Log10FI of reference feature'), ylabel('Log10FI residuals')

else 
    M2S_figureH(0.45,0.35); set(gcf,'Name','Trends')
    plot(refSetToPredict,targetSetToPredict - refSetToPredict,'.k'), axis tight, grid on
    y_lim=ylim; x_lim = xlim; hold on; plot(refSetToPredict,Residuals_trendline-refSetToPredict ,'or','MarkerSize',4), grid on
    ylim(y_lim); xlim(x_lim); xlabel('reference'), ylabel('target - reference')

    M2S_figureH(0.45,0.35); set(gcf,'Name','Residuals')
    plot(refSetToPredict,Residuals_X,'.k'), axis tight, hold on, xlim1 = xlim; grid on
    %plot(xlim',[0;0],'-k')
    xlabel('reference'), ylabel('target - reference') 
end