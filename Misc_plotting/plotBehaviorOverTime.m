function PlotBehaviorOverTime(folders,varargin)
% this function will take a folder or cell of folder and return the animals
% behavioral performance over time for all ARTs found in said folder
%

p

if ~exist('folders','var') || isempty(folders)
    folders = uigetdir;
end


if ~exist('last_n_days','var')
    last_n_days = 'all'

    

if ischar(folders)
    
    folders = {folders};
end

for exp =1:length(folders)
    parent_path = folders{exp};
    matfiles = dir(char(fullfile(parent_path, 'ART*.mat')));
    matfiles = struct2cell(matfiles);
    matfiles = matfiles(1,:)';
    
    %create all unique levels
    uLevels = [];
    h = {};
    for mf_idx = 4:length(matfiles)
        h_temp = WF_getPsignalInfo(char(fullfile(parent_path,matfiles(mf_idx))));
        
        if ~strcmp(strip(h_temp.Class),'Tone')
            continue
        end
        
        
        h{end+1} = h_temp;
        uLevels = unique(cat(1,uLevels,h{end}.uLevels)) ;
    end
    
    m = length(uLevels);
    % loop over said levels and dates to create indicies
    Hits = cell(m,1);
    Miss = cell(m,1);
    Early = cell(m,1);
    Early_idx = cell(m,1);
    Hits_idx = cell(m,1);
    Miss_idx = cell(m,1);
    
    
    trial_idx = 0;
    for mf_idx = 1:length(h)
        
        
        handles = h{mf_idx};
        total_trials = length(handles.Levels);
        [~,a]  =fileparts(handles.Psignalfile);
        fprintf( '%d:  %s , %d Hits,%d Early,%d Trials \n ',mf_idx,...
            a,sum(handles.Hits), sum(handles.Early),total_trials) ;
        
        if sum(handles.Hits) == 0 || total_trials < 15 || length(handles.uLevels)<2 
            continue
        end
        trial_idx = trial_idx+1 ;
            for lvl = 1:length(handles.uLevels)
                lvl_idx = find(handles.uLevels(lvl) == uLevels);
                idx = handles.Levels == handles.uLevels(lvl);
                % values
                Hits_temp  =  handles.Hits(idx);
                hits_idx_temp = ones(1,length(Hits_temp)) *trial_idx;
                pHit(trial_idx,lvl_idx) = sum(Hits_temp)/length(Hits_temp) ;
                
                
                Hits{lvl_idx} =cat(2,Hits{lvl_idx},Hits_temp);
                Hits_idx{lvl_idx} = cat(2,Hits_idx{lvl_idx},hits_idx_temp);
                
                
                Early_temp  =  handles.Early(idx);
                early_idx_temp = ones(1,length(Early_temp)) *trial_idx;
                
                Early{lvl_idx} =cat(2,Early{lvl_idx},Early_temp);
                Early_idx{lvl_idx} = cat(2,Early_idx{lvl_idx},early_idx_temp);
                pEarly(trial_idx,lvl_idx) =  sum(Early_temp)/length(Early_temp);
                
                Miss_temp  =  handles.Miss(idx);
                miss_idx_temp = ones(1,length(Miss_temp)) *trial_idx;
                
                Miss{lvl_idx} =cat(2,Miss{lvl_idx},Early_temp);
                Miss_idx{lvl_idx} = cat(2,Miss_idx{lvl_idx},miss_idx_temp);
                pMiss(trial_idx,lvl_idx) = sum(Miss_temp)/length(Miss_temp) ;
                
            end
        end
        
        
        % package all into behavetime
    if ~exist('pHit','var')
        continue
    end 
    n = size(pHit,1);
    % if we want the last_n_days calculate the indicies of thedays
    % else collect all of the days
    
    if  isnumeric(last_n_days) && last_n_days < n
        
        n_days = n-last_n_days +1 ;
        
    else
        n_days =1;
    end
    % sort for pColor
    pHit = pHit(n_days:end,:)';
    
    
    pEarly =  pEarly(n_days:end,:)';
    
    dPrime = pHit - pEarly;
    
    % Plot
    figure;
    subplot(2,1,1)
    pcolor_fn(pHit,uLevels)
    title(parent_path)
    subplot(2,1,2)
    pcolor_fn(dPrime,uLevels)
    % clear vars
    clear pHit
    clear uLevels
    clear pEarly
    clear pMiss
    
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






function pcolor_fn(P,uLevels)

P = cat(1,P, zeros(1, size(P,2)) );
P = cat(2, P,zeros(size(P,1),1)) ;
pcolor(P)
yticks(.5:1:length(uLevels)+1)
yticklabels(num2str([0;uLevels]))

end













