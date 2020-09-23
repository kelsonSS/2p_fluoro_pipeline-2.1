%CreateFilesForUpload- Script to upload files 
BehaviorFiles = dir('Z:\Kelson\Analyzed\**\2020*\Tones*Active');
PassiveBehaviorFiles =  dir('Z:\Kelson\Analyzed\**\2020*\Tones*Passive');
AllBehaviorFiles = cat(1,BehaviorFiles,PassiveBehaviorFiles);


for file_idx = 1:length(AllBehaviorFiles)
   f_len = fprintf('extracting %s,\n',base_path);
   base_path = fullfile(AllBehaviorFiles(file_idx).folder,...
                        AllBehaviorFiles(file_idx).name);
                    
                  
%repeat for active and passive data
% extractFunction(Passive)
% extractFunction(Active)

% make file Struture
 f_parts = strsep(base_path,filesep);
 f_parts{3} = 'ExtractedBehaviorData';
 out_path = strjoin(f_parts,filesep);
 
 if ~exist(out_path,'dir')
     mkdir(out_path)
 end 
 % Find and copy Psignal File
 try
 psignal_fname = dir([base_path, '\' '*Phys*.mat']);

     
 psignal_fname = psignal_fname.name;
 copyfile(fullfile(base_path,psignal_fname),...
          fullfile(out_path,psignal_fname));
 catch
 end 
      
  % Copy Experiment File
  expt_file = fullfile(base_path,'Fluorescence.mat');
  if exist(expt_file,'file')
      copyfile(expt_file, fullfile(out_path,'Fluorescence.mat'))
  end 
  
  % pupil
  pupil_dir = fullfile(base_path,'pupil'); 
  if exist(pupil_dir,'dir')
      copyfile(pupil_dir,fullfile(base_path,'pupil'))
  end 
 
   f_len = fprintf(repmat('\b',[1,f_len]));
end 

 

                    
                    