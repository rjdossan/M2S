%% Create RTMZ and a label for a dataset
% This creates a label of the type "SLPOS_1250.1234_5.1234" (MS data m/z and RT in minutes)
% or of the type "CPMG_5.1234" (NMR data)
%
% INPUT:
% MS DATA: 
%
% Give a label, and the values for the two columns:
% [X.VarInfo.MZRT_str] = M2S_createLabelMZRT('SLPOS',X.VarInfo.mzmed,X.VarInfo.rtmed)
%
% NMR DATA: 
% 1. Give a label, and the ppm values
% [X.VarInfo.MZRT_str] = M2S_createLabelMZRT('CPMG',X.VarInfo.ppm)
% 
% NOTE1: for string labels, use instead:
% X.VarInfo.MZRT_str = strcat("BILISA","_",stringArray)
%
% NOTE2: THE SECOND STRING ARGUMENT IS DIVIDED BY 60 IN CASE IT IS > 100


function [MZRTstring] = M2S_createLabelMZRT(stringForLabel,varargin)

% If there is only one argument, create only one string
arg1 = varargin{1};
arg1_str=cell(length(arg1),1); 
for a=1:length(arg1)
    
    if isnan(arg1(a))
        arg1_str{a,1} = 'NaN';
    else
        temp_arg1=num2str(floor(10000*arg1(a)));
        if arg1(a)<0.1
            temp_arg1=num2str(floor(100000*arg1(a)));
            arg1_str{a,1} = ['0.0',temp_arg1(end-3:end-1)];
        elseif arg1(a)<1
            arg1_str{a,1} = ['0.',temp_arg1(end-3:end)];
        else
            arg1_str{a,1} = [temp_arg1(1:end-4),'.',temp_arg1(end-3:end)];
        end
    end
end


% if there are two arguments, create a second string    
if length(varargin) == 2   
    arg2 = varargin{2};
    
    %% force RT in minutes
    if nanmax(arg2)>100
        arg2=1/60*arg2;
    end
    
    arg2_str=cell(length(arg2),1); 
    for a=1:length(arg2)
        
        if isnan(arg2(a))
            arg2_str{a,1} = 'NaN';
        else
            temp_arg2=num2str(floor(10000*arg2(a)));
            if arg2(a)<0.0001
                arg2_str{a,1} = '0.0001';
            elseif arg2(a)<0.001      
                temp_arg2=num2str(floor(100000*arg2(a)));
                arg2_str{a,1} = ['0.000',temp_arg2(end:end-1)];
            elseif arg2(a)<0.01      
                temp_arg2=num2str(floor(100000*arg2(a)));
                arg2_str{a,1} = ['0.00',temp_arg2(end-2:end-1)];
            elseif arg2(a)<0.1
                temp_arg2=num2str(floor(100000*arg2(a)));
                arg2_str{a,1} = ['0.0',temp_arg2(end-3:end-1)];
            elseif arg2(a)<1
                arg2_str{a,1} = ['0.',temp_arg2(end-3:end)];
            else
                arg2_str{a,1} = [temp_arg2(1:end-4),'.',temp_arg2(end-3:end)];
            end
        end
    end
end
   
% Create the final label
MZRTstring=cell(length(arg1),1);  
if length(varargin) == 1
    for a=1:length(arg1)
        MZRTstring{a,1}=[stringForLabel,'_', arg1_str{a,1}]; 
    end
elseif length(varargin) == 2
    for a=1:length(arg1)
        MZRTstring{a,1}=[stringForLabel,'_', arg1_str{a,1},'_', arg2_str{a,1}]; 
    end
end

