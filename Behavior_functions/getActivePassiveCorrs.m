
function [active_corr, passive_corr,TNBehavior] = getActivePassiveCorrs(TNBehavior,BehaviorType,Timing)
 
expts = size(TNBehavior,1);
 
if ~exist('Timing','var')
    Timing = 'Tone';
end 

for expt = 1:expts
    
    passive = TNBehavior{expt,1} ;
    active =  TNBehavior{expt,2};
    
   
    excited = squeeze(max(nanmean(active.DFF(30:end,:,:)))) > 0;
    %passive_idx = passive.active{:,2}>0; 
    %active_idx  = active.active{:,2}>0;
    passive_idx = passive.responsive;
    active_idx = active.responsive;
    
    uLevels_expt =sort(unique(active.FreqLevelOrder{:,2}),'descend');
    
    if length(passive_idx) ~= length(active_idx) || sum( active_idx & passive_idx) == 0
         TNBehavior{expt,1} = {} ;
         TNBehavior{expt,2} = {};
         continue
    end 
    responsive_idx = squeeze(passive_idx & active_idx & excited);
    responsive = table([1:length(passive_idx)]',responsive_idx,'VariableNames',{'Neuron','Activity'});
    
    passive.active = responsive;
    active.active = responsive;
   
    switch BehaviorType
        case 'All'
            behave_idx = 1:length(active.handles{1}.Hits);
        case 'Hit'
            behave_idx = find(active.handles{1}.Hits);
          
        case 'Hits'
             behave_idx = find(active.handles{1}.Hits);
            
        case 'Miss'
              behave_idx = find(active.handles{1}.Miss);
        
        case 'Early'
            behave_idx = find(~active.handles{1}.Early);
              
        case 'Incorrect'
            behave_idx = find(~active.handles{1}.Hits);
        otherwise 
            error('Unrecognized behave index')
    end
    behave_idx = behave_idx( behave_idx<= size(active.DFF,2));
    
    active.DFF = active.DFF(:,behave_idx,:);
    active.DFF_Z = active.DFF_Z(:,behave_idx,:);
    active.FreqLevelOrder = active.FreqLevelOrder(behave_idx,:);
    
    if ~isempty( active.FreqLevelOrder)
        active_freq = active.FreqLevelOrder{1,1};
    else
        active_freq = 0;
    end 
    passive_freqs = unique(passive.FreqLevelOrder{:,1});
    
  
     
    [~,passive_idx] =  min(abs(passive_freqs - active_freq));
    PassiveFreq = passive_freqs(passive_idx);
    
    
   
    
    passive_freq_idx = any(passive.FreqLevelOrder{:,1} == PassiveFreq,2);
    passive.DFF = passive.DFF(:,passive_freq_idx,:);
    passive.DFF_Z = passive.DFF_Z(:,passive_freq_idx,:);
    passive.FreqLevelOrder = passive.FreqLevelOrder(passive_freq_idx,:);
    passive.BestFrequencies = BestFrequencyAnalysis(passive);

    
     TNBehavior{expt,1} = passive ;
     TNBehavior{expt,2} = active;
end 

expts_to_use = cellfun(@(x) ~isempty(x),TNBehavior(:,1));

TNBehavior = TNBehavior(expts_to_use,:);

passive_corr = CorrelationsByAnimalBehavior(TNBehavior(:,1),[],Timing);
active_corr = CorrelationsByAnimalBehavior(TNBehavior(:,2),[70,60,50],Timing);
