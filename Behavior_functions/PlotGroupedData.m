function PlotGroupedData(Behavior)

% Takes grouped Data from MungeBehaviorGroupData and plots it 

 Levels = fieldnames(Behavior);
 
 
 HitRate =[];
 LickLatency = [];
 
 figure
 % latency plot 
 m = length(Levels);
 for curr_level = 1:m
       subplot(m,1,curr_level);
       histogram(Behavior.(Levels{curr_level}).SNRLatency,...
           'BinWidth',.1);
       hold on 
       title(Levels{curr_level},'Interpreter','None')
       ax = axis;
       
       plot([1,1],[ax(3), ax(4)], 'g-') 
       plot([2,2],[ax(3), ax(4)], 'g-') 
              
    
     
 end 
 
 
 % SNR Plot 1 -   HitRate grouped by SNR
 
 hitRate_mu = zeros(m,1);
 hitRate_SEM = zeros(m,1);
 hitRateAll= [];
 
 for curr_level = 1:m
      hitRate_mu(curr_level) = nanmean(Behavior.(Levels{curr_level}).HitRateMean);
      hitRate_SEM(curr_level) = nanstd(Behavior.(Levels{curr_level}).HitRateMean) /...
                     sqrt(length(Behavior.(Levels{curr_level}).HitRateMean)) ;
      try
        hitRateAll = cat(2,hitRateAll,Behavior.(Levels{curr_level}).HitRateMean); 
      catch
          temp = zeros(length(hitRateAll),1);
          hit_lvl= Behavior.(Levels{curr_level}).HitRateMean;
          for ii = 1:length(hit_lvl)
              temp(ii) = hit_lvl(ii);
          end 
           hitRateAll = cat(2,hitRateAll,temp(1:length(hitRateAll)));
           clear temp
      end 
          
 
 end             
     figure
 bar(hitRate_mu)
 set(gca,'TickLabelInterpreter','None')
 hold on
 xticks(1:m);
 xticklabels(Levels);
 errorbar(hitRate_mu,hitRate_SEM ,'.')
 title('Fraction Correct vs. dB SNR')
 ylabel('Fraction Correct')
 xlabel ('dB SNR')
       

               
 
 % SNR Plot 2 -HitRate  grouped by Animal 
 levelLabels = strrep(Levels,'SNR_','');
 levelLabels = strrep(levelLabels,'minus_','-');
 figure

 hold on
 shadedErrorBar([],mean(hitRateAll),std(hitRateAll)/sqrt(6),'k')
  plot(hitRateAll')
 xticks(1:4)
 set(gca,'Xdir','rev')
 xticklabels(levelLabels)
 ylabel('Percent Hit Rate')
xlabel('dB SNR') 
 
 