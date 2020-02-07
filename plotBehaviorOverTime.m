function PlotBehaviorOverTime(folders)
% this function will take a folder or cell of folder and return the animals 
% behavioral performance over time for all ARTs found in said folder 
% 
if ~exist('folders','var')
    folders = uigetdir;
end 

if ischar(folders)
    
    folders = {folders};
end

for exp = 1:length(folders)
    parent_path = folders{exp};
    matfiles = dir(char(fullfile(parent_path, 'ART*.mat')));
    matfiles = struct2cell(matfiles);
    matfiles = matfiles(1,:)';
    
    %create all unique levels
    uLevels = [];
    for mf_idx = 1:length(matfiles)
       h{mf_idx} = WF_getPsignalInfo(char(fullfile(parent_path,matfiles(mf_idx))));
       uLevels = unique(cat(1,uLevels,h{mf_idx}.uLevels)) ;  
    end
    
    % loop over said levels and dates to create indicies 
    
    for mf_idx = 1:length(matfiles)
        handles = h{mf_idx};
       
        if length(handles.uLevels)>1 
            for lvl = 1:length(handles.uLevels)
                idx = handles.Levels == handles.uLevels(lvl);
                % values
                Hits_temp  =  handles.Hits(idx);
                Early_temp  =  handles.Early(idx);
                Miss_temp  =  handles.Miss(idx);
                % indicies
                hits_idx_temp = ones(1,length(Hits_temp)) *mf_idx;
                early_idx_temp = ones(1,length(Early_temp)) *mf_idx;
                miss_idx_temp = ones(1,length(Miss_temp)) *mf_idx;
                
                Hits{lvl} =cat(2,Hits{lvl},Hits_temp);
                Early{lvl} =cat(2,Early{lvl},Early_temp);
                Miss{lvl} =cat(2,Miss{lvl},Early_temp);
                
                Hits_idx = cat(2,Hits_idx{lvl},hits_idx_temp);
                Early_idx = cat(2,Early_idx{lvl},early_idx_temp);
                Miss_idx = cat(2,Miss_idx{lvl},miss_idx_temp);
                
            end
        end 
                      
                         % package all into behavetime
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                    
                    end 
        
%         
%         for ul = 1:uLevels
%             test_level =  uLevels(ul);
%             lvl_idx = find(handles.uLevels == test_level);
%             if lvl_idx
%                 % extract hits and early trials (misses can be found
%                 % from non-hits, non-early trials as there are no FAs)
%                 hits_temp  = handles.HitsLevels{lvl_idx};
%                 hits_temp_idx = ones(1,length(hits_temp));
%                 early_temp = handles.Early{lvl_idx};
%                 early_temp_idx = ones(1,length(early_temp));
%                 
%                 BehaveTime.Hits{test_level} = cat(2,BehaveTime.HitsLevels{test_level},...
%                     hits_temp);
%                 BehaveTime.Early{test_level} = cat(2,BehaveTime.EarlyLevels{test_level},...
%                     hits_temp);
%                 BehaveTime.HitsIdx{test_level} = cat(2,BehaveTime.HitsIdx{test_level},hits_temp_idx);
%                 BehaveTime.EarlyIdx{test_level} = cat(2,BehaveTime.EarlyIdx{test_level},early_temp_idx);
%             end
            
            
            
            
end
      
        




   % plotting 









end

    


    
    
    
    
    
    

