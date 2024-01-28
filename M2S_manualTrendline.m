% [mousePointCoords, predictedPchip] = M2S_manualTrendline(x_pointsToPredict)
% This function is part of M2S toolbox to match metabolomics features across untargeted datasets.
% 
% Function to manually define x,y points, create a curve using them and
% predict the y-coordinates of a set of new x-values.
%
% INPUT:
% x_pointsToPredict: x-points vector to predict. Can be used empty as: M2S_manualTrendline()
%
% OUTPUT:
% mousePointCoords: the coordinates chosen with the mouse
% predictedPchip: y-coordinates predicted using polynomial cubic interpolation
%
% NOTE: the function acts in any of the plots of the current figure (gcf)
%
%   Rui Climaco Pinto, 2022
%   Imperial College London
%   Adapted from https://uk.mathworks.com/matlabcentral/answers/128400-draw-points-using-mouse
%   Thanks Jon!

function [mousePointCoords, predictedPchip] = M2S_manualTrendline(x_pointsToPredict)



fprintf('\nFunction to define a trendline in a plot.\n')
fprintf('The last figure to be clicked is the currently active figure (gcf).\n')
fprintf('Choose multiple points in the plot, starting from the left, using mouse left-click.\n')
fprintf('To finish click ENTER.\n')
fprintf('To cancel click Ctrl+C\n')

fig = ancestor(gca,'figure');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate array to hold mouse click coordinates
maxPoints=2000;
mousePointCoords = zeros(maxPoints,2);
% set up loop to collect and display mouse click points
count = 0;
for k = 1:maxPoints
    % get the mouse click point, or terminate if user presses enter
    %  in which case the coordinates will be returned empty
    coords = ginput(1);
    if isempty(coords)
        break
    end
    count = count + 1;
    mousePointCoords(count,:) = coords;
    
    plot(mousePointCoords(k,1),mousePointCoords(k,2),'or','MarkerSize',8); hold on
    if k > 1
        plot(mousePointCoords(k-1:k,1),mousePointCoords(k-1:k,2),'-r','MarkerSize',8);        
    end
end
% clean up
hold off
mousePointCoords = mousePointCoords(1:count,:); % trim off unused array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Predict points using polynomial cubic interpolation. Plot

if nargin == 0
    predictedPchip = pchip(mousePointCoords(:,1),mousePointCoords(:,2),mousePointCoords(:,1));
    hold on, 
    plot(mousePointCoords(:,1),predictedPchip,'*b')
elseif nargin == 1
    predictedPchip = pchip(mousePointCoords(:,1),mousePointCoords(:,2),x_pointsToPredict);
    hold on, 
    plot(x_pointsToPredict(:,1),predictedPchip,'*b')
end














































