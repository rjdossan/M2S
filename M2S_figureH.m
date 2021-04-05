%% Create or position a figure in the centre of the screen
% To create:
% figureH = M2S_figureH(widthPercent,heightPercent)
% To position
% figureH = M2S_figureH(widthPercent,heightPercent,figH)
%
% widthPercent and heightPercent are between 0 and 1

function figureH = M2S_figureH(widthPercent,heightPercent,figH)

%Ensure root units are pixels and get the size of the screen:
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

%Define the size and location of the figures:
pos1  = [round(0.5*(1-widthPercent) * scnsize(3)),... 
    round(0.5*(1-heightPercent) *scnsize(4)),...
    round(widthPercent * scnsize(3)),...
    round(heightPercent * scnsize(4))];

if nargin==2
    %Create the figures:
    figure('OuterPosition',pos1) 
    figureH=gcf;
    
else
    figureH = figH;
    set(figH,'OuterPosition',pos1) 
end
