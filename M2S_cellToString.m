%% Xstring = M2S_cellToString(Xcell,tempChar)
% Function to transform a cell variable to a string variable
%
% It works for variables in a table, and when there are empty matrices in cells
% e.d. in LCMS annotations e.g., {[]}
% The default function call substitutes empty cells for '.' and later
% deletes the '.'. If another char is preferred to the '.' (e.g. because
% there are positions in the cell with '.') then add it as 'tempChar'

function Xstring = M2S_cellToString(Xcell,tempChar)
if nargin == 1
    tempChar = '.';
end

if iscell(Xcell)
    Xcell(cellfun(@isempty,Xcell)) = {tempChar};
    Xstring = string(Xcell);
    Xstring(Xstring==tempChar) = "";
else
    Xstring = Xcell;
end



