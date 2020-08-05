function idx = FluorescenceFileCheck(RawFile,in_path,save_path)

% this function takes a RawFile location along with the base path and save
% path and checks if a corresponding Fluorescence file was made, indicating
% that this file has been extracted. Part of TwoPhotonPipeline 
% KS 2020 


   
   in_path_length = length(in_path);
   
   RawFile = RawFile(in_path_length+1:end);
   
   file_to_test  = fullfile(save_path,RawFile,'Fluorescence.mat');
  %
   idx = exist(file_to_test,'file');
   
   

    