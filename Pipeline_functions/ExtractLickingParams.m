function ExtractLickingParams(input,debug)




if ~exist('debug','var')
    debug= 0;
end 

    TimingInfo = struct('Errors',[]);
    %% path prep
 
    PsignalFileName = input.psignalfiles;
    SavePath = input.path;
    LocalPath = input.path;
        
    TimingFile = fullfile(LocalPath,'TimingInfo.mat');
    xmlfile = fullfile(LocalPath,'Experiment.xml');
    PsignalFile= fullfile(LocalPath,PsignalFileName);
      fps =input.expectedFPS ;
     if exist(char(fullfile(LocalPath,'redchannel.mat')),'file')
         fps = fps/2;
     end 
     
     try
     xml = get_options_from_xml(xmlfile);
     catch
       warning('bad xml file')
         return
     end 
     
    %% File Prep
    if iscell(PsignalFile)
        PsignalFile = PsignalFile{1};
    end 
    
     if ~exist(PsignalFile,'file')
        TimingInfo.Errors = LocalPath;
         return
    end
    
    
   % if file was already processed correctly, continue to next expt 
    if exist(TimingFile,'file') 
            load(TimingFile)
    else 
                return
    end 
    
    ThorSyncFile = fullfile(LocalPath, 'Episode001.h5');
    h5 = h5info(ThorSyncFile);
    
    if ~contains([h5.Group(2).Datasets.Name],'Licks')
        return
    end 
    
    
    
    licks_raw=h5read(ThorSyncFile,'/AI/Licks');
    
    licks = findpeaks(licks,5.1)
    
    licks = unique(floor(licks/10))* 10; 
    
    
    

    

        

