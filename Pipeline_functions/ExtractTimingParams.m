function TimingInfo =  ExtractTimingParams(input,debug)

if ~exist('debug','var')
    debug= 0;
end 

animalID = input.animalID;
for e = 1:length(input.expname)
    
    %% path prep
    if ~isfield(input,'manual_mode') || input.manual_mode ~= 1
    ExpName = input.expname{e};
    PsignalPath = input.psignalfiles{e};
    SavePath = fullfile(input.savepath,animalID);
    LocalPath = fullfile(input.path,ExpName);
    else 
        ExpPath = input.expname{1};
       [SavePath, ExpName] = fileparts(ExpPath);
        LocalPath = ExpPath;
      PsignalPath =  dir(fullfile(ExpPath,'RND*.mat'));
      if isempty(PsignalPath)
          PsignalPath =  dir(fullfile(ExpPath,'ART*.mat'));
      end
      
      PsignalPath = PsignalPath.name;
    
    end
    
    TimingFile = fullfile(SavePath,ExpName,'TimingInfo.mat');
    xmlfile = fullfile(LocalPath,'Experiment.xml');
    PsignalFile=fullfile(SavePath , ExpName, PsignalPath);
      fps =input.expectedFPS ;
     if exist(fullfile(SavePath,ExpName,'redchannel.mat'),'file')
         fps = fps/2;
     end 
     
     xml = get_options_from_xml(xmlfile);
     
    %% File Prep
    try
        if ~exist( fullfile(SavePath,ExpName,PsignalPath),'file' )
            copyfile(fullfile(LocalPath, ExpName, PsignalPath ),...
                     fullfile(SavePath, ExpName, PsignalPath) );
        end
        
     
    catch
        
     
        continue
    end
    
    
   
    if exist(TimingFile,'file') && debug ~= 1 
            load(TimingFile)
            continue 
    end 
    %% get actual number of frames in registered file 
    
    fh = fopen( fullfile(SavePath,ExpName,'greenchannelregistered.raw'));
    if fh == -1
        continue
    end 
    
   num_frames_actual = FindRawImgSize(fh,[xml.dimX xml.dimY]);
   fclose(fh);
        
        
        %% extraction    
        
        
        ThorFile = fullfile(LocalPath, 'Episode001.h5');
        if exist(ThorFile,'file')
            % ThorImage >3.1
            if  ~exist(fullfile(SavePath,ExpName, 'Episode001.h5'),'file')
                copyfile( ThorFile , fullfile(SavePath,ExpName) )
                
            end
            
            TimingInfo = getTimingInfo_H5(ThorFile, PsignalFile, fps,xml,num_frames_actual);
     
            
        else
            
            ThorFile = fullfile(LocalPath ,'timing.txt' );
            % Thorimage 2.1 
            if exist(ThorFile,'file')
                if ~exist(fullfile(SavePath,ExpName, 'timing.txt'),'file')
                    
                    copyfile( ThorFile , fullfile(SavePath,ExpName) )
                end
                
                TimingInfo = getTimingInfo(ThorFile, PsignalFile, fps);
               
            end
        end


if  isfield(TimingInfo,'FrameIdx') % if extraction was successful

if debug
    % if there is already a file check that these changes dont affect the
    % data
    if exist(TimingFile,'file')
        old  = load(TimingFile);
        old = old.TimingInfo;
        % if there is a shift more than 1 frame stop for manual inspection
      try
        shift_flag = max(max(abs(old.FrameIdx - TimingInfo.FrameIdx))) > 1; 
      catch
          shift_flag = 1;
      end 
        
        if shift_flag == 1 
             fo=h5read(ThorFile,'/DI/Frame Out');
             try
                 gate=h5read(ThorFile,'/AI/ai5'); % trial gate signal
             catch
                 gate=h5read(ThorFile,'/AI/PsignalGate');
             end
            fprintf('Trial 1 New:%d-%d \n Trial 1 Old:%d-%d \n',...
                   TimingInfo.FrameIdx(1,:),old.FrameIdx(1,:))
             % plot frame out and trial gate 
               figure; plot(gate); hold on; plot(fo);
             
             title(sprintf('%s \n %s',  SavePath,input.expname{e}),...
                 'interpreter','none')
               xlabel(sprintf('Trial 1 New:%d-%d \n Trial 1 Old:%d-%d',...
                   TimingInfo.FrameIdx(1,:),old.FrameIdx(1,:) ));
             clear gate 
             clear fo
        end 
    end 
end 

% save the timing file 

save(TimingFile,'TimingInfo')
clear TimingInfo
end
end 

end 