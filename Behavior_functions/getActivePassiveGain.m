function Out = getActivePassiveGain(TNBehavior,SaveName)

if ~exist('SaveName','var')
    SaveName = []
end 
   
    passive_mu_all = [];
    passive_norm_all = [];
    active_mu_all = [];
    active_norm_all = [];
    Gain_levels = struct('SNR_20dB',[],'SNR_10dB',[],'SNR_0dB',[]);
    

for expt = 1:length(TNBehavior)

    active = TNBehavior{expt,2};
    %passive = subsetPassive(TNBehavior{expt,1},active);
    passive = TNBehavior{expt,1};
    
    
if   size(passive.DFF,3) ~= size(active.DFF,3)
    fprintf('expt %d mismatch \n',expt)
    continue
end 

    
    r_idx = passive.responsive  & active.responsive;
    
    clean_idx = passive.Clean_idx & active.Clean_idx;
    
    final_idx = r_idx & clean_idx;
    
    DF_mu_passive = squeeze(nanmean(passive.DFF(:,:,final_idx),2));
    DF_mu_active = squeeze(nanmean(active.DFF(:,:,final_idx),2));
    DF_norm_passive = squeeze(nanmean(passive.DFF_norm(:,:,final_idx),2));
    DF_norm_active = squeeze(nanmean(active.DFF_norm(:,:,final_idx),2));
    [Gain_levels] = getActiveLevels(active,DF_mu_passive,final_idx,Gain_levels);

   % seperate into active-Hit and active-Miss and active-early ? 
    
    passive_mu_all = cat(2,passive_mu_all,DF_mu_passive);
    passive_norm_all = cat(2,passive_norm_all,DF_norm_passive);
    active_mu_all = cat(2,active_mu_all,DF_mu_active);
    active_norm_all = cat(2,active_norm_all,DF_norm_active);
    
    
        
end 

 Gain = active_mu_all(1:90,:) - passive_mu_all;

   neg_idx = max(Gain(40:70,:)) <= 0;
%  neg_idx = max(passive_mu_all(40:70,:)) <= 0;
%  neg_idx = max(active_mu_all(40:70,:)) <= 0;
  sum(neg_idx)
  length(neg_idx)
  figure

  subplot(1,3,2)
  active_order=PlotMeanResponse(active_mu_all(1:90,:));

  title('Active')
 subplot(1,3,1)
  PlotMeanResponse(passive_mu_all,active_order);
  title('Passive')
  subplot(1,3,3)
  PlotMeanResponse(Gain,active_order);
  title('Gain')
 
  saveas(gcf,[SaveName '-Gain.pdf'])
  
  
figure 
plotShadedErrorBar(Gain(:,neg_idx),'b')
hold on 
plotShadedErrorBar(Gain(:,~neg_idx),'r')
ylim([-15 5])


title( 'Attentional Gain')
xlabel( 'Time(s)')
xticks(0:30:90)
xticklabels(0:1:3)
ylabel('DFF (active - passive)') 
 saveas(gcf,[SaveName '-AverageGain.pdf'])



Out.NegCell_Gain = max(Gain(40:70,neg_idx))';
Out.PosCell_Gain = max(Gain(40:70,~neg_idx))';
% seperate in to + and n
PlotCI(Out.NegCell_Gain)
PlotCI(Out.PosCell_Gain,gcf)
legend({'Neg-Gain','','Pos-Gain'})
ylim([-15 5])

if SaveName
    saveas(gcf,[SaveName '-AttentionalGainBar.pdf'])
end

Out.GainLevels = PlotLevelGain(Gain_levels);
if SaveName
    saveas(gcf,[SaveName '-AttentionalGain_Levels.pdf'])
end









function GainLevels = getActiveLevels(active,passive,final_idx,GainLevels)

 Levels = [70,60,50];
 LevelNames = {'SNR_20dB','SNR_10dB','SNR_0dB'};
 
 for lvl = 1:length(Levels)
     
     lvl_idx = active.FreqLevelOrder{:,2} == Levels(lvl);
     if any(lvl_idx)
        DF_mu_active = squeeze(nanmean(active.DFF(:,lvl_idx,final_idx),2));
        Gain_lvl = DF_mu_active(1:90,:) - passive;
        GainLevels.(LevelNames{lvl}) = cat(2,GainLevels.(LevelNames{lvl}), Gain_lvl); 
     end
 end
 
 
 
 
 function gain_levels = PlotLevelGain(gain_levels)
        
    
f = fieldnames(gain_levels);
n_levels = length(f);
figure 

for field_idx = 1 : n_levels
    fname = f{field_idx};    
    
    Gain = gain_levels.(fname);  
    
    
    neg_idx = max(Gain(40:70,:)) <= 0;
    
    subplot(1,n_levels,field_idx)
    plotShadedErrorBar(Gain(:,neg_idx),'y')
    hold on 
    plotShadedErrorBar(Gain(:,~neg_idx),'y')
    ylim([-20 10])
    xticks(0:30:90)
    xticklabels(1:3)
    xlabel('Time (S)')
    title(fname,'Interpreter','none')
    
    
end 

for field_idx = 1 : n_levels
    
    fname = f{field_idx};    
    
    Gain = gain_levels.(fname);  
    
    
    neg_idx = max(Gain(40:70,:)) <= 0;
    
    pos_gain = mean(Gain(40:70,~neg_idx));
    neg_gain  = mean(Gain(40:70,neg_idx));
    
    
    gain_levels.pos_gain_levels{field_idx} = pos_gain;
    gain_levels.neg_gain_levels{field_idx} = neg_gain;
    
      
    
end 






 
 
 





function Cell_order = PlotMeanResponse(DFF2,Cell_order)
            

%% sort by average  max response 

[~, DFF_order] = sort(max(DFF2));
DFF2 = DFF2(:,DFF_order);
DFF2 = DFF2(:,:);
%% plot average responses 
%imagesc(DFF2')
%FramesToSeconds
%caxis([-20 20])
%% Sort by time to peak and normalize responses to the max response 
[DFF2_max, DFF2_timing]  = max(DFF2);
% normalize
DFF3 = DFF2./DFF2_max; 


% sort
if ~exist('Cell_order','var')
[DFF3_timing,Cell_order]=  sort(DFF2_timing);
end 
DFF3 = DFF3(:,Cell_order);

% plot
imagesc(DFF3')
caxis([-1 1 ])
FramesToSeconds
title('Time To Max Response')




function FramesToSeconds 
ax = gcf;
figure(ax) 
xticks([30,60,90,120])
xticklabels({'1','2','3','4'})
xlabel('Time (sec)')
ylabel('Neurons')



function subsetPassive(Passive,active)


detect_freq = unique(active.handles{1}.FreqLevelOrder{:,1});

passive_all_freqs = unique(passive.handles{1}.FreqLevelOrder{:,1});
passive_freq_idx = find( min(abs(passive_Freqs -...
                         detect_freq)))
   
                     



