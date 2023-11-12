%% Function to apply to a table to transform all cell variables into string 
% NOTE: First load the data using the import tool in matlab or "readtable"
% Then run the function:
% tableFinal = M2S_tableCellToString(tableInitial)
%
% Rui Pinto 2021

function tableFinal = M2S_tableCellToString(tableInitial)
tableInitial_classes = varfun(@class,tableInitial,'OutputFormat','cell');
cells_idx = find(string(tableInitial_classes) == 'cell');
tableFinal = table; 
for c=1:size(tableInitial,2)
    if sum(cells_idx == c)>0
        tableFinal = [tableFinal, table(cellToString(table2array(tableInitial(:,c))),'VariableNames',tableInitial.Properties.VariableNames(c))];
    else
        tableFinal = [tableFinal, table(table2array(tableInitial(:,c)),'VariableNames',tableInitial.Properties.VariableNames(c))];
    end
end
tableFinal_classes = varfun(@class,tableFinal,'OutputFormat','cell');
%disp(array2table([string(tableInitial_classes'),string(tableFinal_classes')],'VariableNames',{'initialClasses','finalClasses'}));

   