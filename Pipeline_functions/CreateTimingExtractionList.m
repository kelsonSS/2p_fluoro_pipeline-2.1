function [NeedsExtraction,PsignalFiles] = CreateTimingExtractionList(savepath)
% this function looks to see what files are ready to be proccessed and
% extracted via the TwoPhotonPipeline's ExtractFluorescence functionality

if ~exist('savepath','var')
    savepath = 'Z:\kelson\analyzed';
end 
 
 Registered = dir( [savepath, '\**\greenchannelregistered.raw']);
 Registered = struct2cell(Registered);
 Registered = Registered(2,:)';
 
NeedsExtraction_idx = false(length(Registered),1);
for  file_idx = 1:length(Registered) 
  if isempty( dir([Registered{file_idx}, '\**\TimingInfo.mat']) )
       NeedsExtraction_idx(file_idx) = 1;
  end
end
  
 NeedsExtraction = Registered(NeedsExtraction_idx);
 
 
 PsignalFiles = cellfun(@PsignalFileCheck,NeedsExtraction,'UniformOutput',0);

 
end 


