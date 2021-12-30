function Out = getFluoroHitsMiss(TNBehavior,out_folder)
    
if ~exist('out_folder','var')
    to_plot = 0
else
    to_plot = 1
    
    if ~exist(out_folder,'dir')
        mkdir(out_folder)
    end 
end 



noise_level = 50;
detectLevelSNR = [ 0 ,10,20];

Out = CreateOutputStructure(detectLevelSNR);


detectLevels = detectLevelSNR  + noise_level;    

    for expt_idx =1:length(TNBehavior)
        
        
        expt= TNBehavior{expt_idx};
        handles = expt.handles{1};
        uLevels = handles.uLevels - 50;
       
        
        F = expt.DFF(:,:,expt.Clean_idx & expt.active{:,2}>0);
        %F = expt.DFF_Z(:,:,expt.Clean_idx & expt.active{:,2}>0);
        
        hits_idx = handles.Hits';
        miss_idx = handles.Miss';
        early_idx = handles.Early';
        Levels = handles.Levels;
       
        all_lvl_idx = any(handles.Levels - noise_level == [ 0 10 20],2);
        
        
        n_expts = size(F,2)
         if length(all_lvl_idx) ~= n_expts
             continue
         end 
%            hits_idx = hits_idx(1:n_expts);
%            miss_idx = miss_idx(1:n_expts);
%            early_idx = early_idx(1:n_expts);
%            all_lvl_idx = all_lvl_idx(1:n_expts);
%            Levels = Levels(1:n_expts)
%         end
        
        
        all_hits_temp = squeeze(nanmean(F(:,all_lvl_idx & hits_idx,:),2));
        all_miss_temp = squeeze(nanmean(F(:,all_lvl_idx & miss_idx,:),2));
        all_early_temp = squeeze(nanmean(F(:,all_lvl_idx & early_idx,:),2));
        
        
        
        Out.All.Hits = cat(2,Out.All.Hits,all_hits_temp);
        Out.All.Miss = cat(2,Out.All.Miss,all_miss_temp);
        Out.All.Early = cat(2,Out.All.Early,all_early_temp);

        
        for lvl = 1:length(uLevels)
            uLevel = uLevels(lvl);
            if any(uLevel == detectLevelSNR)
                
                lvl_idx = (Levels - noise_level) == uLevel;
                hits_temp = squeeze(nanmean(F(:,lvl_idx & hits_idx,:),2));
                miss_temp = squeeze(nanmean(F(:,lvl_idx & miss_idx,:),2));
                early_temp = squeeze(nanmean(F(:,lvl_idx & early_idx,:),2));
             
                
                
                
                snr_name = sprintf('SNR_%d',uLevel);
                snr_name = strrep(snr_name,'-','minus_');
               
                Out.(snr_name).Hits =cat(2,Out.(snr_name).Hits,hits_temp);
                Out.(snr_name).Miss =cat(2,Out.(snr_name).Miss,miss_temp);
                Out.(snr_name).Early = cat(2,Out.(snr_name).Early,early_temp);

            end 
        end 
    end 
    
    if to_plot
        PlotResults(Out,out_folder)
    end 

    
end 



function Out =CreateOutputStructure(levels)

levels = [num2cell(levels),'All'];

for ii = 1:length(levels)
if isnumeric(levels{ii})
    snr_name = sprintf('SNR_%d',levels{ii});
    snr_name = strrep(snr_name,'-','minus_');
else 
    snr_name = levels{ii};
end 

Out.(snr_name).Hits = [];
Out.(snr_name).Miss = [];
Out.(snr_name).Early = [];
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
