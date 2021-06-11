function [NeedsExtraction,PsignalFiles] = CreateTimingExtractionList(savepath)
% this function looks to see what files are ready to be proccessed and
% extracted via the TwoPhotonPipeline's ExtractFluorescence functionality

if ~exist('savepath','var')
    savepath = '\\vault3\data\kelson\analyzed';
end 
 
 Registered = dir( [savepath, '\**\greenchannelregistered.raw']);
 Registered = struct2cell(Registered);
 Registered = Registered(2,:)';
 
 TimingInfo = dir( [savepath, '\**\TimingInfo.mat']);
 TimingInfo = struct2cell(TimingInfo);
 TimingInfo = TimingInfo(2,:)';

 

 NeedsExtraction =  setdiff(Registered,TimingInfo);

 
 
 PsignalFiles = cellfun(@PsignalFileCheck,NeedsExtraction,'UniformOutput',0);

 
end 


