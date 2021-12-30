function getBehavioralPerformanceFromHandles(Experiment_files) 

% this function grabs the experimental performance from one or more psignal
% files given to it. this function will then show both individual
% performance and group average performance for said files 

% Kelson Shilling-Scrivo 2019 
 


for fl_idx  = 1:length(Experiment_files) 

   handles = WF_getPsignalInfo(Experiment_files{fl_idx}) ;
   
   uF = unique(handles.Freqs);
   uL = sort(unique(handles.Levels),'descend'); 
   m =  length(uL);
  %% init 
  hits = cell(m,1);
  miss = cell(m,1);
  early = cell(m,1);
  for lvl = 1:m 
   %% gather level data 
   
   expt_idx =  handles.Levels == uL(lvl) ;    
   
   hits{lvl} =  handles.Hits(expt_idx);
   miss{lvl} = handles.Miss(expt_idx);
   early{lvl} = handles.Early(expt_idx); 
   
   hitrate(lvl) = sum(hits{lvl} == 1) / length(hits{lvl}) ; 
   missrate(lvl) = sum(miss{lvl} == 1) / length(miss{lvl}) ; 
   earlyrate(lvl) = sum(early{lvl} == 1) / length(early{lvl}) ; 
   
  end
  %% save for group data 
 hits_all{1}(:,fl_idx) = hits;
 miss_all{1}(:,fl_idx) = miss ;
 early_all{1}(:,fl_idx) = early;

 hitrate_all(:,fl_idx) = hitrate;
 missrate_all(:,fl_idx) = missrate;
 earlyrate_all(:,fl_idx) = earlyrate;
 


end 
%% plot individual data




end 
