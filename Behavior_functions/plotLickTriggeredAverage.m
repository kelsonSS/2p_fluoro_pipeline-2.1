function plotLickTriggeredAverage(TNBehavior)

if ~exist('pad_size','var')
    pad_size = 5;
end 

EarlyLick_DFFs = []; 
HitLick_DFFs = [];
    for expt_idx =1:length(TNBehavior)
        
        
        expt= TNBehavior{expt_idx};
        
        
        if isempty(expt.LickTiming)
            continue
        end
        
        data = getDataExpt(expt,EarlyLick_DFFs,pad_size);
        
        
    end  
    
    figure
    plotShadedErrorBar(data,'k')
    shifted_xticks = arrayfun(@num2str,[1:size(data,1)]-pad_size,'UniformOutput',0)
    hold on 
    xlim([1 10])
    xticks(5)
    xticklabels('0')
    ylabel('DFF')
    xlabel('time rel. lick')
    plot([5 5],[0 20],'k--')
    title("Lick Triggered Average")
    
end

function data = getDataExpt(e,data,pad_size)
  trial_data = [];
  handles = e.handles{1};
  early_idx = logical(handles.Early);
  lick_times = e.LickTiming.TrialLickFrames;
  %get early trials  
  if ~any(early_idx)
      return
  end 
    
    fr_frames =cellfun(@getFirstResponse,lick_times);
    fr_frames = fr_frames(early_idx);
    
    early_fluoro = e.DFF(:,early_idx,:);
    
    
    
    trial_data = getLickWindow(early_fluoro,fr_frames,pad_size);
    
    data = cat(2,data,trial_data);
    
  
  
    
  
  
end 


function d = getLickWindow(DFF,Licks,pad_size)

% pad array and adjust Licks accordingly
DFF = padarray( DFF , [0, pad_size, 0], 'pre');

Licks = round(Licks + pad_size);

DFF = DFF( Licks-pad_size:Licks + pad_size,:,:);

d = reshape(DFF,size(DFF,1),[]);



end 

function r = getFirstResponse(x)

if isempty(x)
    r = nan;
else 
    r = x(1);
end 

end 


function PlotResults(Out,out_path)

f = fieldnames(Out);

for lvl =1:length(f)
    
 pos_idx = max(Out.(f{lvl}).Hits(40:70,:)) > 0
 hits_pos = Out.(f{lvl}).Hits(:,pos_idx);
 hits_neg = Out.(f{lvl}).Hits(:,~pos_idx);
 miss_pos = Out.(f{lvl}).Miss(:,pos_idx);
 miss_neg = Out.(f{lvl}).Miss(:,~pos_idx);
 early_pos = Out.(f{lvl}).Early(:,pos_idx);
 early_neg = Out.(f{lvl}).Early(:,~pos_idx);
 figure; hold on
shadedErrorBar([],nanmean(hits_pos,2),nanstd(hits_pos,[],2)/ sqrt(length(hits_pos)) * 1.96, 'b')
shadedErrorBar([],nanmean(miss_pos,2),nanstd(miss_pos,[],2)/ sqrt(length(miss_pos)) * 1.96 )
%shadedErrorBar([],nanmean(early_pos,2),nanstd(early_pos,[],2)/ sqrt(length(early_pos)) * 1.96,'r')
title_name = sprintf('%s-pos',f{lvl});
title(title_name,'Interpreter','none');
%ylim([-1 5])
print(fullfile(out_path,[title_name '.pdf']),'-dpdf')

 figure; hold on
shadedErrorBar([],nanmean(hits_neg,2),nanstd(hits_neg,[],2)/ sqrt(length(hits_neg)) * 1.96, 'b')
shadedErrorBar([],nanmean(miss_neg,2),nanstd(miss_neg,[],2)/ sqrt(length(miss_neg)) * 1.96 )
%shadedErrorBar([],nanmean(early_neg,2),nanstd(early_neg,[],2)/ sqrt(length(early_neg)) * 1.96,'r')
title_name = sprintf('%s-neg',f{lvl});
title(title_name,'Interpreter','none');
%ylim([-5 1])
print(fullfile(out_path,[title_name '.pdf']),'-dpdf')





 figure
 shadedErrorBar([],nanmean(hits_pos - miss_pos,2),...
                   nanstd(hits_pos - miss_pos,[],2)/ sqrt(length(hits_pos)) * 1.96)
 title_name =sprintf('%s: Hit-miss-pos',f{lvl});
% ylim([-5 5])
  title(title_name ,'Interpreter','none');
 
 
 print(fullfile(out_path, sprintf('%s_HitMissDiff-pos.pdf',f{lvl}) ),'-dpdf' ) 

  figure
 shadedErrorBar([],nanmean(hits_neg - miss_neg,2),...
                   nanstd(hits_neg - miss_neg,[],2)/ sqrt(length(hits_neg)) * 1.96)
 title_name =sprintf('%s: Hit-miss-neg',f{lvl});
  title(title_name ,'Interpreter','none');
 % ylim([-5 5])
 
 print(fullfile(out_path, sprintf('%s_HitMissDiff-neg.pdf',f{lvl}) ),'-dpdf' ) 


 

end 
end 





