function Out = getActivePassiveGain(TNBehavior,SaveName)

if ~exist('SaveName','var')
    SaveName = []
end 
   
    passive_mu_all = [];
    passive_norm_all = [];
    active_mu_all = [];
    active_norm_all = [];


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
 
figure 
plotShadedErrorBar(Gain(:,neg_idx),'b')
hold on 
plotShadedErrorBar(Gain(:,~neg_idx),'r')
plotShadedErrorBar(Gain(:,:),'k')
title( 'Attentional Gain')
xlabel( 'Time(s)')
xticks(0:30:90)
xticklabels(0:1:3)
ylabel('DFF (active - passive)') 

Out.NegCell_Gain = Gain(70,neg_idx)';
Out.PosCell_Gain = Gain(70,~neg_idx)';
% seperate in to + and n
PlotCI(Gain(70,neg_idx)')
PlotCI(Gain(70,~neg_idx)',gcf)
legend({'Neg-Gain','','Pos-Gain'})
ylim([-15 5])

if SaveName
    saveas(gcf,[SaveName '-AttentionalGainBar.pdf'])
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
   
                     



