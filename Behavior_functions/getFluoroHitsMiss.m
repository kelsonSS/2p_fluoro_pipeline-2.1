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
        
        all_hits_temp = correctBaseline(all_hits_temp,5);
        all_miss_temp = correctBaseline(all_miss_temp,5);
        all_early_temp = correctBaseline(all_early_temp,5);
        
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


function  DF_Corrected = correctBaseline(DF,baseline_frames)

% baseline correction
trialdur = size(DF,1);

B_Vec = repmat(nanmean(DF(1:baseline_frames,:)),[trialdur,1]);
DF_Corrected = DF -B_Vec;




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
 
 lvl_name = strrep(f{lvl},'-','_');
 % plot pos
 title_name = sprintf('%s_pos',lvl_name);
 PlotFigs(hits_pos,miss_pos,early_pos,title_name,[-2 18])
 PlotFigs(hits_pos(1:30,:),miss_pos(1:30,:),early_pos(1:30,:),['prestim_' title_name],[-3 4])
% plot neg
 title_name = sprintf('%s_neg',lvl_name);
 PlotFigs(hits_neg,miss_neg,early_neg,title_name,[-20 5])
 PlotFigs(hits_neg(1:30,:),miss_neg(1:30,:),early_neg(1:30,:),['prestim_' title_name],[-8 6]) 


%% old figs 
% 
%  figure
%  shadedErrorBar([],nanmean(hits_pos - miss_pos,2),...
%                    nanstd(hits_pos - miss_pos,[],2)/ sqrt(length(hits_pos)) * 1.96)
%  title_name =sprintf('%s: Hit-miss-pos',f{lvl});
% % ylim([-5 5])
%   title(title_name ,'Interpreter','none');
%  
%  
%  print(fullfile(out_path, sprintf('%s_HitMissDiff-pos.pdf',f{lvl}) ),'-dpdf' ) 
% 
%   figure
%  shadedErrorBar([],nanmean(hits_neg - miss_neg,2),...
%                    nanstd(hits_neg - miss_neg,[],2)/ sqrt(length(hits_neg)) * 1.96)
%  title_name =sprintf('%s: Hit-miss-neg',f{lvl});
%   title(title_name ,'Interpreter','none');
%  % ylim([-5 5])
%  
%  print(fullfile(out_path, sprintf('%s_HitMissDiff-neg.pdf',f{lvl}) ),'-dpdf' ) 
% 

 

end 

PlotDiffFig(Out,f)

end 


function PlotDiffFig(Out,lvls)

hits_high = Out.(lvls{1}).Hits;


end 


function PlotFigs(hits,miss,early,title_name,y_lim)
if ~exist('y_lim','var')
    y_lim = [];
end 
hits = correctBaseline(hits,10);
miss = correctBaseline(miss,10);
early = correctBaseline(early,10);


ShadedErrorBarFig(hits,miss,early,title_name,y_lim)
BarPlotFig(hits,miss,early,title_name)

end 

function ShadedErrorBarFig(hits,miss,early,title_name,y_lim)

figure; hold on
 

shadedErrorBar([],nanmean(hits,2),nanstd(hits,[],2)/ sqrt(length(hits)) * 1.96, 'g')
shadedErrorBar([],nanmean(miss,2),nanstd(miss,[],2)/ sqrt(length(miss)) * 1.96 ,'k' )
shadedErrorBar([],nanmean(early,2),nanstd(early,[],2)/ sqrt(length(early)) * 1.96,'m')
title(title_name,'Interpreter','none');
if y_lim
    ylim(y_lim)
end 
%ylim([-1 5])
xticks(0:30:120)
xticklabels(0:4)
xlabel('Time (S)')
   
print([title_name '.pdf'],'-dpdf')
end 

function BarPlotFig(hits,miss,early,title_name)

end 
