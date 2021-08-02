function out = ComparePassiveActive(TNBehavior)


Difference_means ={};
for expt = 1:size(TNBehavior,1)
    
    
    Passive = TNBehavior{expt,1};
    Active = TNBehavior{expt,2};
    
    active_Freq = Active.handles{1}.uFreqs;
    Passive_mean = Mean_DFF(Passive,active_Freq);
    Active_mean = Mean_DFF(Active,active_Freq);
    % ensure sizes are compatible 
    Active_mean = Active_mean(1:size(Passive_mean,1),:);
    
    Difference_mean = Active_mean - Passive_mean;
    
    Difference_means{expt} = Difference_mean;
    
end 

  
  figure
    plotShadedErrorBar(cell2mat(Difference_means))



end 

function DFF_mu = Mean_DFF(experiment,uFreqs)

handles = experiment.handles{1};
if ~exist('uFreqs','var')
    uFreqs = handles.uFreqs;
end

trials = find(any(handles.Freqs == uFreqs',2));

trials = trials( trials <= size(experiment.DFF,2));

DFF_mu = squeeze(nanmean(experiment.DFF_norm(:,trials,:),2));



end 