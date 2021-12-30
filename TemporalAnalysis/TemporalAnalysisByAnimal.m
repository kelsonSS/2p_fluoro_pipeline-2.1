function [Timings]= TemporalAnalysisByAnimal(DF,lvl,SaveName)

if ~exist('SaveName','var')
    SaveName = ''
end 

bin_size = 10;

Cell_stats = TemporalAnalysis(DF,lvl,SaveName);

 % create indicies  
if exist('lvl','var') && ~strcmp('all',lvl)
    lvl_idx = any(DF.FreqLevelOrder{:,2} == lvl , 2);
    assert(sum(lvl_idx) >0 , 'Incorrect level used') 
else
    lvl_idx = true(1,size(DF.DFF,2));
end 
    

 
  
   active_idx=  DF.active{:,2}>0  ;
     
   % get average reponse
   DFF = squeeze(nanmean(DF.DFF(1:end-2,lvl_idx,active_idx),2));
    
    % get bins 
   hist_bins = [1:bin_size:size(DFF,1)+4];


   % get expt_idx 
    expt_list = DF.experiment_list(active_idx);
    
    
    n_expts = max(expt_list);
    
    
Tone_on_off_all = [];
Noise_on_off_all = [];
    
    timing_col_id = 1;
    for expt_id = 1:n_expts
      st = GetSoundTimesFromHandles(DF.handles{expt_id});
 
       DFF_temp = DFF(:, expt_list == expt_id);
   
        if size(DFF_temp,2) < 10 
          continue
        end 
       % get first derivative 
   
   
   DFF_temp = diff(DFF_temp);
   
   [dff_max,timing] =  max(DFF_temp);
%      Tone_on_max =   max(DFF_temp(st.Tone_on:st.Tone_on+bin_size,:));
%      Tone_off_max =  max(DFF_temp(st.Tone_off:st.Tone_off+bin_size,:));
%      Noise_on_max =  max(DFF_temp(st.Noise_on:st.Noise_on+bin_size,:));
%      Noise_off_max = max(DFF_temp(st.Noise_off:st.Noise_off+bin_size,:));
%  
%      Tone_on_off_ratio = Tone_on_max - Tone_off_max)/ (Tone_on_max + Tone_off_max;
%      Noise_on_off_ratio = Noise_on_max ./ Noise_off_max;
%      
%      Tone_on_off_all = cat(2,Tone_on_off_all,Tone_on_off_ratio);
%      Tone_on_off_expts{expt_id} = Tone_on_off_ratio;
%  
%      Noise_on_off_all = cat(2,Noise_on_off_all,Noise_on_off_ratio);
%      Noise_on_off_expts{expt_id} = Noise_on_off_ratio;
% 



    



    exc_idx = dff_max > 0;
    timing = timing(exc_idx);
   

    timing_prc(timing_col_id,:) = histcounts(timing,hist_bins)./length(timing);

    timing_col_id = timing_col_id +1;
    end 
    
   % Tone on/offset
    ON_time = find(hist_bins == st.Tone_on);
    OFF_time =find(hist_bins == st.Tone_off);
  % noise onset
    Noise_ON_time = find(hist_bins == st.Noise_on);
    %Tone_ratio_animal = timing_prc(:,ON_time)./(timing_prc(:,ON_time) + timing_prc(:,OFF_time));
    Tone_bias_animal = ( timing_prc(:,ON_time) - timing_prc(:,OFF_time) ) ./...
                       ( timing_prc(:,ON_time) + timing_prc(:,OFF_time) );
    Tone_balance_number = sum( abs(Tone_bias_animal) < .5)  ;       
    Tone_balance_percentage  =  Tone_balance_number /...
                              length(Tone_bias_animal) * 100;            

    Noise_bias_animal = ( timing_prc(:,Noise_ON_time) - timing_prc(:,ON_time) ) ./...
                        ( timing_prc(:,Noise_ON_time) + timing_prc(:,ON_time) );


    try
   [~,main_effects,stats] = anova1(timing_prc, [] ,'off');
  timing_stats= multcompare(stats,'Display','off');
  timing_stats=timing_stats(1:size(timing_prc,2),:)  
    catch
    end 


Timings.main_effects = main_effects;
Timings.Stats =  timing_stats;
Timings.Stats_Cell = Cell_stats;
Timings.TimingByAnimal_prc = timing_prc;
Timings.ToneBiasByAnimal = Tone_bias_animal;
Timings.NoiseBiasByAnimal = Noise_bias_animal;
Timings.ToneBalanceNumber = Tone_balance_number;
Timings.ToneBalancePercentage = Tone_balance_percentage;
%Timings.ToneRatiosByAnimal = Tone_ratio_animal;
% Timings.ToneRatios_expt = Tone_on_off_expts;
% Timings.NoiseRatios_all = Noise_on_off_all;
% Timings.NoiseRatios_expt = Noise_on_off_expts;


    %% plotting
   
    
    
    figure; bar(mean(timing_prc),'BarWidth',1)
    hold on 
    errorbar(mean(timing_prc), std(timing_prc) / sqrt(n_expts) * 1.96 ,'.')
    
    % plotting individual animals with scatter plot
    
    %make a matrix denoting each timebin 
    
    n_timebins = size(timing_prc,2);
    n_expts = size(timing_prc,1);
    
    timing = repmat( [1:n_timebins],n_expts,1);
    
 scatter(timing(:),timing_prc(:),15, 'MarkerEdgeColor','k','MarkerFaceColor','w')
 ylim([0 .7])
  ticks = yticks;
 yticklabels(ticks * 100)
 ylabel('%')
title([SaveName '-Timing'],'Interpreter','None')
 if SaveName
     saveas(gcf,sprintf('%s-TimingByAnimal.pdf', SaveName))
 end
 
figure 
histogram(Tone_bias_animal,[-1:.25:1])
ax = gca
hold on 
plot([ -.5 , -.5],[ax.YLim(1),ax.YLim(2)],'k--')
plot([  .5 , .5],[ax.YLim(1),ax.YLim(2)],'k--')
title([SaveName '-ToneBias'],'Interpreter','None' )
xticks([-1:.25:1])
xlabel('ToneBias')
ylabel('Count')
if SaveName
     saveas(gcf,sprintf('%s-ToneBiasAnimal.pdf', SaveName))
 end




    

           
         
