function [NeedsExtraction,PsignalFiles] = CreateCellExtractionList(savepath)
% this function looks to see what files are ready to be proccessed and
% extracted via the TwoPhotonPipeline's ExtractFluorescence functionality


if ~exist('savepath','var')
    savepath = '\\vault3\data\kelson\analyzed';
end 

TimingInfo = dir( [savepath, '\**\TimingInfo.mat']);
 TimingInfo = struct2cell(TimingInfo);
 TimingInfo = TimingInfo(2,:)';
 
 CellDefinitions = dir( [savepath, '\**\CellDefinitions.mat']);
 CellDefinitions = struct2cell(CellDefinitions);
 CellDefinitions = CellDefinitions(2,:)';
 
 Fluorescence = dir( [savepath, '\**\Fluorescence.mat']);
 Fluorescence = struct2cell(Fluorescence);
 Fluorescence = Fluorescence(2,:)';
 
 
 HasCellDef = intersect(TimingInfo,CellDefinitions);
 NeedsExtraction =  setdiff(HasCellDef,Fluorescence);

 
 
 PsignalFiles = cellfun(@PsignalFileCheck,NeedsExtraction,'UniformOutput',0);

 
end 


function  psignal_files = PsignalFileCheck(Directory) 


PathFiles =  dir(Directory);
PathFiles = struct2cell(PathFiles);
PathFiles = PathFiles(1,:);
% find files
expression = '(ART)|(RND)'; 
matches =  regexp( PathFiles,expression,'once');

% convert logicnl index to row number
psignal_file_idx = find(~cellfun(@isempty,matches));

% extract file 
psignal_files = PathFiles{psignal_file_idx};



end 


