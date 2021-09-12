function plotDFTracesLevels(Data,activity)
% this function takes Sorted Neural Traces in the form: 
%   [Frame X Trial]
% and creates the average response for each condition

% Data denotes Data in the Format produced by Fluoro_to_Table
% activity is the activity as recorded from active_list 

% note- for FRAs Expts should be sorted by frequency first and then by level 
% requires ToneInNoise_MeanTrace which is a stylized shadedErrorBar

if ~exist('activity','var')
    activity = 3;
end
% init

FreqLevels = sortrows(unique(Data.FreqLevelOrder),{'Freqs','Levels'},{'Ascend','Descend'});
T = length(unique(Data.FreqLevelOrder{:,1})); 
L = length(unique(Data.FreqLevelOrder{:,2})); 
num_subplots = T * L;
active_list = Data.active{:,2};
active_idx = active_list >= activity;
cell_ids = find(active_idx);

DFTraces = Data.DFF(:,:,active_idx);
 

[~,Max_idx]=  sort(squeeze(max(mean(abs(DFTraces),2))),'descend');
DFTraces = DFTraces(:,:,Max_idx);
cell_ids = cell_ids(Max_idx);

Trials = size(Data.DFF,2);
TrialsPerFreq = floor( Trials / ( T * L ) ) ;







% main loop
for cell = 1:size(DFTraces,3) 
figure('Position', [20 50 1800 900])
Start = 1;
FL_idx = 1; 
suptitle(sprintf('Neuron: %d ', cell_ids(cell)) )
for ii = 1:T
    for jj = 1:L
         idx = ii + (T)*(jj-1); % plots in column major order
       % find traces 
       
        Freq  = FreqLevels{FL_idx,1};
        Level = FreqLevels{FL_idx,2};
       trial_idx = Data.FreqLevelOrder{:,1} == Freq & ...
                   Data.FreqLevelOrder{:,2} == Level ;
               
       
       
       
       % plot trace
        subplot(L,T,idx)
         hold on
        
        
        ToneInNoise_MeanTrace(DFTraces(:,trial_idx,cell));
        title( sprintf( 'Freq: %d, level: %d',Freq, Level))
       
        %move to next condition
        FL_idx = FL_idx+1; 
        
          axis off
end 


end  
ax = findall(gcf,'Type','Axes')
ax = ax(1:end-1);
ymin = min([ax.YLim]);
ymax = max([ax.YLim]);

set(ax,'Ylim',[ymin ymax])

% style subplots 
for idx  = 1:num_subplots
    subplot(L,T,idx)
style_subplot()
 end 


pause
 end 
  

function style_subplot()
 aa = axis;
  plot([30 30], [aa(3),aa(4)],'--')
  plot([60 60], [aa(3),aa(4)],'--g')
  plot([90 90], [aa(3),aa(4)],'--g')
  plot([120 120], [aa(3),aa(4)],'--r')
  xlim([0 150])       
        
        
        


