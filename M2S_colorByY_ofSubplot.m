
function M2S_colorByY_ofSubplot(subplotNr,figH)
load ('M2ScolorScheme.mat');

if nargin == 1
    figH = gcf;
end

% Find number of subplots
n=numel(figH.Children);

% Find how many subplot rows and columns
axx = findobj(figH,'type','axes');
pos = cell2mat(get(axx,'position'));
nrows = numel(unique(pos(:,2))); % the same Y position means the same row
ncols = numel(unique(pos(:,1))); % the same X position means the same column

% Get the y values of the desired subplot, to use as colours
subplot(nrows,ncols,subplotNr)
ax = gca; 
h = findobj(gca,'Type','line');
x = (h(1).XData)'; 
y = (h(1).YData)'; 

hold on, scatter(x,y,20*ones(size(x)),y,'filled')
colormap(M2Scolormap)

pause(1)
%% Colour the other subplots using the same colours as above
not_subplotNr = setdiff(1:n,subplotNr);

for s = 1:length(not_subplotNr);
    subplot(nrows,ncols,not_subplotNr(s))
    ax = gca; 
    h1 = findobj(gca,'Type','line');
    x1 = (h1(1).XData)'; 
    y1 = (h1(1).YData)'; 
    hold on, scatter(x1,y1,20*ones(size(x1)),y,'filled')
end
