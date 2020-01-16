function b = get_behavioral_performance(file)


%%
%Get list of animals interested
if ~exist('file','var')
    [fl dir] = uigetfile('MultiSelect','on');
    
    file = fullfile(dir,fl);
end
% get last N trials from that animal

%%

% For each animal extract their hit, miss and False alarm statistics
% by level
if ~iscell(file)
    temp =  {file};
    file = temp;
    clear temp
end
   
handles = WF_getPsignalInfo(file{1});
m = length(handles.uLevels);
    levels_prev = {};
    hits = zeros(1,m); 
    misses = zeros(1,m); 
    false_alarms =zeros(1,m) ;
    n_trials = zeros(1,m); 



for ii = 1:length(file)
    
    handles= WF_getPsignalInfo(file{ii}) ;
    
    disp(ii)
    %% test that all levels are similar to previous expts12
    levels = handles.uLevels;
    
%     if  ii> 1  && ~isempty(setxor(levels_prev,levels))
%         warning('these files do not have the same levels')
%     end
    levels_prev = levels;
    
    %% extract level behavioral performance
    m = length(handles.uLevels);
    n = length(handles.Hits);
    
    if n < 10
        continue
    end 
    for uL = 1:m
        level_idx = handles.Levels == handles.uLevels(uL);
        
        hits(uL) = hits(uL) + sum(handles.Hits(level_idx));
        misses(uL) = misses(uL) + sum(handles.Miss(level_idx));
        false_alarms(uL) = false_alarms(uL) + sum(handles.Early(level_idx));
        n_trials_temp(uL) = n_trials(uL) + sum(level_idx);
     
     
        
    end
    
 
 n_trials = n_trials + n_trials_temp;
 fprintf('%d,%d,%d',ii, sum(n_trials),length(handles.Hits))
end


%% packaging
b.hits = hits;
b.miss = misses;
b.false_alarm = false_alarms;
b.n_trials = n_trials
b.levels = handles.uLevels - handles.OverallDB;
b.dir = file;
