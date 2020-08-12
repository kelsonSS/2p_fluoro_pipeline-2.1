function[NeedsExtraction,PsignalFiles,animalIDs] = createDataList(input)

% this function checks the folder containing every raw file in 
% \\VAULT3\Data\Kelson and checks to see if there is a corresponding 
% folder in ..\Kelson\Analyzed. if there is not that indicates that the
% file needs to be analyed and creates the appropriate Datalist.
%

% munging 
if isfield(input,'inpath')
    file_path = input.inpath;
else     
    file_path = '\\VAULT3\Data\Kelson\Files to unload' ;
end 

if isfield(input,'savepath')
    save_path = input.savepath ;
else 
    save_path = '\\vault3\data\kelson\analyzed';
end 

image_name = input.regexp;


% querying
 RawImages = dir( [file_path, '\**\', image_name ]);
 RawImages = struct2cell(RawImages);
 RawImages = RawImages(2,:)';
 
 PsignalFiles = cellfun(@PsignalFileCheck,RawImages,'UniformOutput',0);
 % Find 
 HasPsignalFiles = ~cellfun(@isempty,PsignalFiles);
 HasFluorescence = cellfun(@(x)... 
                    FluorescenceFileCheck(x,file_path,save_path),...
                    RawImages);
 
 RawImages = RawImages(HasPsignalFiles & ~HasFluorescence );
 PsignalFiles = PsignalFiles(HasPsignalFiles & ~HasFluorescence );
 
 
 NeedsExtraction = cellfun(@(x) RelativePathFromFullPath(x,file_path),...
     RawImages,'UniformOutput',0);
 animalIDs = cellfun(@(x) AnimalIDFromPath(x,file_path),...
     RawImages,'UniformOutput',0);
 
 
end 


    
    
   
    
    
 
     

