 function filePreparation(expt_path,load_path,save_path,psignalfile_name)
        % this is legacy code designed to create a new path for the extracted data
        % and move the appropriate
        % legacy code from Nik Francis 2016
        
      
            % To-do start fixing below
            %Local path for data
            old_path = fullfile(load_path,expt_path) ;
            new_path = fullfile(save_path,expt_path) ;
            
            if ~exist(new_path,'file')
                mkdir(new_path)
            end 
            
            % move pupil Data
            
            if exist(fullfile(old_path,'pupil'),'dir')
                     copyfile(fullfile(old_path,'pupil'),...
                     fullfile(new_path,'pupil'))
            end 
            % Copy psignalfile 
            
            copyfile(fullfile(old_path,psignalfile_name),...
                     fullfile(new_path,psignalfile_name));
            
            
            %Move relevant ThorImage files to local path
            Thor_path= fullfile(old_path,...
                'Experiment.xml');
            
            copyfile(Thor_path,fullfile(new_path,'Experiment.xml'));
 
            % there should be either a thor timing file or H5 file depending on
            % using 2.1 or 3.0 respestively we move either one over here
            try
                Timing_path =fullfile(old_path,...
                    'timing.txt');
                if ~exist( fullfile(new_path, 'timing.txt'),'file')
                    copyfile(Timing_path, fullfile(new_path,'timing.txt'));
                end
            catch
              H5 =   dir(fullfile(old_path, '**/*/Episode001.h5') );
                if ~isempty(H5)
                    h5file_path = fullfile(H5.folder,H5.name);
                    try
                    movefile(h5file_path,fullfile(new_path, 'Episode001.h5'));
                    catch
                    warning('H5 File Not moved')    
                    end 
                end
            end
        end
        
  

