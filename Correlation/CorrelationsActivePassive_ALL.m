function out = CorrelationsActivePassive_ALL(TNBehavior,SaveName,PlotColor,Timing)

if ~exist('PlotColor','var') 
    PlotColor = 'k'
end 

if ~exist('Timing','var')
    Timing = 'Tone'
end 

if ~exist('SaveName','var')
    SaveName = [];
end 


[active_All,passive_All,TNBehavior] = getActivePassiveCorrs(TNBehavior,'All',Timing);
active_Hits = getActivePassiveCorrs(TNBehavior,'Hit',Timing);
active_Miss = getActivePassiveCorrs(TNBehavior,'Miss',Timing);
active_Early = getActivePassiveCorrs(TNBehavior,'Early',Timing);
active_Incorrect = getActivePassiveCorrs(TNBehavior,'Incorrect',Timing);



close_idx_all  = CompareBestFrequencyToBehaviorTarget(TNBehavior);

close_idx = cell2mat(close_idx_all);


out.Hits =  FullCompareCorrs(passive_All,active_Hits,...
                             close_idx,PlotColor,...
                           [SaveName '_Hits']);
                                            


out.All = FullCompareCorrs(passive_All,active_All,...
                          close_idx,PlotColor,...
                          [SaveName '_All'],...
                                out.Hits);


%out.Miss =  FullCompareCorrs(passive_All,active_Miss,...
%                             close_idx,PlotColor,...
%                           [SaveName '_Miss'],...
%                               out.Hits);
                           
out.Early =  FullCompareCorrs(passive_All,active_Early,...
                             close_idx,PlotColor,...
                           [SaveName '_Early'],...
                               out.Hits);

out.Incorrect =  FullCompareCorrs(passive_All,active_Incorrect,...
                             close_idx,PlotColor,...
                           [SaveName '_Incorrect'],...
                               out.Hits);                           

                           
Behaviors = {'All','Hits','Early','Incorrect'};


for b_idx = 1:length(Behaviors)
        
    curr_B = Behaviors{b_idx};
    
    
    out.(curr_B).Corr_by_distance = CorrsByDistance(TNBehavior,10,...
                                                  [SaveName '_' curr_B '_Distance'],...
                                                  curr_B);
                           
end                       



    

        
        
        
        

