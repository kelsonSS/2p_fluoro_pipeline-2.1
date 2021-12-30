function Data =  ReplaceInf(Data,rep_value)

% replaces Inf values in the dataset with rep_value, which is nan by
% default

if ~exist('rep_value','var')
    rep_value = nan;
end 

Data(isinf(Data)) = rep_value;