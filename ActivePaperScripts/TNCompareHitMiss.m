function out =  TNCompareHitMiss(TNBehavior)


out  ={};
Reward_means = {};
Activity_means = {};
Hit_means = {};
Miss_means ={};
FA_means = {};
for expt = 1:size(TNBehavior,1)
    
    Active = TNBehavior{expt,2};
    
    handles = Active.handles{1};
   
    % Gather Behavioral indices
    Hits_idx = handles.Hits;
    Miss_idx = handles.Miss;
    FA_idx   = handles.Early;
    
    
    Hit_mean = GetMeanDFF(Active,Hits_idx);
    Miss_mean = GetMeanDFF(Active,Miss_idx);
    FA_mean = GetMeanDFF(Active,FA_idx);
    
    
     Hit_means{expt} = Hit_mean;
     Miss_means{expt} = Miss_mean; 
    %FA_means = 
   
   
    Reward_mean_expt = Hit_mean - FA_mean;
    Activity_mean_expt = Hit_mean - Miss_mean;
    
   
    Reward_means{expt} = Reward_mean_expt;
    Activity_means{expt} =Activity_mean_expt;
end 

figure
hold on 
 plotShadedErrorBar(cell2mat(Hit_means),'b')
 plotShadedErrorBar(cell2mat(Miss_means),'k')
 StyleFigure()

 figure
    plotShadedErrorBar(cell2mat(Reward_means))
    title('Hit - False Alarm')
    StyleFigure()
 figure
    plotShadedErrorBar(cell2mat(Activity_means))
    title('Hit - Miss')
    StyleFigure()

end 

function DFF_mu = GetMeanDFF(experiment,trials_idx)

active = experiment.active{:,2} > 0;
trials_idx = logical(trials_idx);

trials_idx = trials_idx( 1:size(experiment.DFF,2));

DFF_mu = squeeze(nanmean(experiment.DFF_norm(:,trials_idx,active),2));


end

function StyleFigure()
aa = axis;
hold on
plot([aa(1) aa(2)], [0 0 ] ,'k--')
plot([30 30], [aa(3) aa(4)],'b--')
plot([60 60], [aa(3) aa(4)],'b--')
xlim([0 90])
xticks(0:30:91)
xticklabels({'0','1','2','3'})
xlabel('Time(s)')
ylabel('Normalized Fluorescence(AU)')
    end 



