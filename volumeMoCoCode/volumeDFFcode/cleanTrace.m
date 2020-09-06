%% Algorithmic approach to get rid of erratic 2p traces originally used by Dan Winkowski
% Just made it into an easily applicable function -- Zac B
function [cleanedTrace,keptRois] = cleanTrace(f)
%input trace should be ROIs x Time

roiMin = min(f,[],2);
frameMin = min(f,[],1);
keptRois = (roiMin > 2*median(frameMin));
cleanedTrace = f(keptRois,:);

end