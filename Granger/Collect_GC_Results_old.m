function GC =  Collect_GC_Results(starting_dir) 
% This function takes a collection of output directories from 
% GCAnalysis functions and collects them by condition 
GC.mode = {'hit','miss'};
% TonesLevels-{'TR- Tone', 'TR +20db SNR','TR +10db SNR','TR +0db SNR'};
% hitmiss -  {'Hit','Miss'}
%hitLevels -{'hit-0db','hit-10db','hit-20db'}
  index =  1; 
  
  Run ='Yes' ;
  
 while  strcmp(Run, 'Yes');
    d =  uigetdir(starting_dir);
    d = dir(d);
    d = d(startsWith({d.name},'Data'));
    
    for ii = 1:length(d)
        load( fullfile(d(ii).folder,d(ii).name) )
    end 
    
    % collection
    GC.Exptlist{index} = d(ii).folder;
    GC.GCnumbers(index,:) = GCnumbers';
    GC.GClengths{index,:} = GClengths';
    GC.GCangles{index,:} = GCangles';
    GC.GCStrengthHigh{index} =  GCJH;
    GC.GCStrengthLow {index} =  GCJL;
    
    

    
    
 
Run = questdlg('Continue?');
index = index + 1 ;
end 


