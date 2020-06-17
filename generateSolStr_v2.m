function [solStr,axisStr] = generateSolStr_v2(sol,axis,sensor)
%GENERATESOLSTR returns formatted strings for use in titles
%   generateSolStr(sol,axis) returns a string 'solStr' with the sol(s)
%   covered by the data in 'sol', and an optional axis string 'axisStr' corresponding
%   to the data axis used, specified by 'axis'. 'axis' may be omitted if
%   using a Vector Spectrogram, for example.

if ~exist('axis','var')
    axis = [];
end

if ~exist('sensor','var')
    sensor = [];
end

uniqueSols = unique(sol);
if length(uniqueSols) == 1
    solStr = ['Sol ' num2str(uniqueSols(1))];
elseif length(uniqueSols) == 2
    solStr = ['Sols ' num2str(uniqueSols(1)) '-' num2str(uniqueSols(2))];
end

if ~isempty(axis)
    if axis == 1
        axisStr = 'Z';
    elseif axis == 2
        axisStr = 'N';
    elseif axis == 3
        axisStr = 'E';
    else
        disp('Invalid Axis');
    end
    
    axisStr = [axisStr '\_' sensor];
    
else
    axisStr = [];
end

end

