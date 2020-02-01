function TimingInfo =  ExtractTimingParams(input,debug)

if ~exist('debug','var')
    debug= 0;
end 
%

for e = 1:length(input.expname)
    
    %% path prep
    ExpName = input.expname{e};
    PsignalFile = input.psignalfiles{e};
   
   
    
    [~,animalID] = fileparts(input.path);
    
    SavePath = fullfile(input.savepath,animalID);
    LocalPath = fullfile(input.path,ExpName);
    TimingFile = fullfile(SavePath,ExpName,'TimingInfo.mat');
    xmlfile = fullfile(LocalPath,'Experiment.xml');
      fps =input.expectedFPS ;
     if exist(fullfile(SavePath,ExpName,'redchannel.mat'),'file')
         fps = fps/2;
     end 
     
     xml = get_options_from_xml(xmlfile);
     
    %% File Prep
    try
        if ~exist( fullfile(SavePath,PsignalFile),'file' )
            copyfile(fullfile(input.path, ExpName, PsignalFile ),...
                     fullfile(SavePath, ExpName, PsignalFile) );
        end
        
     
    catch
        
        warndlg('could not move %s %s to new Analyzed folder,Skipping expt',...
             input.path,ExpName)
        warning('could not move %s %s to new Analyzed folder,Skipping expt',...
           input.path,ExpName )
        continue
    end
    
    
    PsignalFile=fullfile(SavePath , ExpName, PsignalFile);
    if exist(TimingFile,'file') && debug ~= 1 
            load(TimingFile)
            continue 
    end 
    %% get actual number of frames in registered file 
    
    fh = fopen( fullfile(SavePath,ExpName,'greenchannelregistered.raw'));
   num_frames_actual = FindRawImgSize(fh,[xml.dimX xml.dimY]);
   fclose(fh);
        
        
        %% extraction    
        
        
        ThorFile = fullfile(LocalPath, 'Episode001.h5');
        if exist(ThorFile,'file')
            
            if  ~exist(fullfile(SavePath,ExpName, 'Episode001.h5'),'file')
                copyfile( ThorFile , fullfile(SavePath,ExpName) )
                
            end
            
            TimingInfo = getTimingInfo_H5(ThorFile, PsignalFile, fps,xml,num_frames_actual);
            save(TimingFile,'TimingInfo')
            
        else
            
            ThorFile = fullfile(LocalPath ,'timing.txt' );
            
            if exist(ThorFile,'file')
                if ~exist(fullfile(SavePath,ExpName, 'timing.txt'),'file')
                    
                    copyfile( ThorFile , fullfile(SavePath,ExpName) )
                end
                
                TimingInfo = getTimingInfo(ThorFile, PsignalFile, fps);
                save(TimingFile,'TimingInfo')
            end
        end


 end 

end 