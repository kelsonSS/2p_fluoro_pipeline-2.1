% activeAnalysis Script 


% preprocessing 
for expt = 1:length(Active)
    root_fl = Active{expt}.DataDirs;
    matfiles = dir([root_fl, '\ART*']);
    matfiles = [matfiles;dir([root_fl,'\RND*'])];
    Active{expt}.PsignalFiles = matfiles.name;
end  % DFF is Time X Trial x Cell

DFF =  Active{1}.DFF;
neuron_max =  max(max( DFF));

% repmat to original size to normalize DFF
max_matrix = repmat(neuron_max,size(DFF,1),size(DFF,2),1);
%size(max_matrix)
% normalized DFF
Active{expt}.DFF_norm = DFF ./ max_matrix; 


Hits = Active{expt}.handles.Hits;
Miss = Active{expt}.handles.Miss;
FA   = Active{expt}.handles.Early;


%% create a plot of a neurons activity be behavioral type 
behavior_modes = {Hits,Miss,FA;'Hits','Miss','False_Alarms'};
n_neurons = size(Active{expt}.DFF,3);


for neuron = 1:n_neurons
     
      % activity selection - plot only active neurons   
      if Active{expt}.active{neuron,2}== 0
          continue
      end 
     
     % figure creation
        figure
  
        hold on 
        
    
    
    for mode = 1:length(behavior_modes)
        
        subplot(1,3,mode)
        
        % selection
        behave_idx =  logical(behavior_modes{1,mode});        
        activity = Active{expt}.DFF_norm(:,behave_idx,neuron); 
        activity_all{mode} = activity;
        % plotting 
        ToneInNoise_MeanTrace(activity,[],'ToneNoiseActive');
              title( sprintf( 'Neuron %d-%s', neuron,behavior_modes{2,mode} ))  
           
    end
    
    % overwrite FA's and put hits - miss average (hack) 
        cla
        hit_miss =   squeeze(mean(activity_all{1},2))...
                   - squeeze(mean(activity_all{2},2));
          plot(hit_miss);
               linkaxes
         
      
    pause
end 





 clear activity
 clear activity_all       

% Clustering 
 k = 3;
for mode = 1:length(behavior_modes)
    
    active_idx = Active{expt}.active{:,2} > 0;
    
    
    % selection
    behave_idx =  logical(behavior_modes{1,mode});
    activity = Active{expt}.DFF_norm(:,behave_idx,active_idx);
    
    % clustering
   
   cluster_idx{mode} =  Cluster_DF(activity,k);
   
end

figure
for cluster = 1:k
    
    
    c_idx =  find(cluster_idx{1} == cluster);
    on_mu = Active{expt}.DFF_norm(:,behave_idx,neuron)
        neuron = c_idx(nn););
    for nn = 1:length(c_idx)
        
        % activity selection - plot only active neurons
        %if Active{expt}.active{neuron,2}< 2
         %   continue
      % end
        
        for mode = 1:2
            subplot(1,2,mode)
            hold on
            % selection
            behave_idx =  logical(behavior_modes{1,mode});
            activity = Active{expt}.DFF_norm(:,behave_idx,neuron);
            ToneInNoise_MeanTrace(activity,[],'ToneNoiseActive')
            title(sprintf('cluster %d- %s',k,behavior_modes{2,mode}))
            linkaxes
        end
        
        
    end
    figure
    % plotting hits v misses
end




% Behavior by level
for expt = 1:length(Active)
  
    uLevels = Active{expt}.handles.uLevels;

    for ul = 1:length(uLevels)
        level_idx = Active{expt}.handles.Levels == uLevels(ul);
        behavior_mat(expt,1,ul) = sum(Active{expt}.handles.Hits(level_idx));
        behavior_mat(expt,2,ul) = sum(Active{expt}.handles.Miss(level_idx));
        behavior_mat(expt,3,ul) = sum(Active{expt}.handles.Early(level_idx));
        
    end
end 


total_trials = sum(behavior_mat);
%% test
for expt = 1:length(Active)
    GCanalTrialsModBalanced_TN( { Active{expt}.DataDirs },'HitMiss')
end
%%
for expt = 1:length(Active)
    GCanalTrialsModBalanced_TN( { Active{expt}.DataDirs },'HitLevels')
end




