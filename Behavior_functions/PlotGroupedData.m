function PlotGroupedData(Behavior,AnimalInfo)

% Takes grouped Data from MungeBehaviorGroupData and plots it 

 Levels = Behavior.uSNRs;
 
 
 HitRate =[];

 figure
 % latency plot 
 m = length(Levels);
 for curr_level = 1:m
       subplot(m,1,curr_level);
        
       histogram(vertcat(Behavior.LickLatency{curr_level,:}),...
           'BinWidth',.1);
       hold on 
       title(Levels(curr_level),'Interpreter','None')
       ax = axis;
       
       plot([1,1],[ax(3), ax(4)], 'g-') 
       plot([2,2],[ax(3), ax(4)], 'g-') 
       xlim([0 , 3])
       xticks(0:1:3)
      
       ylabel('Count')
 end 
 xlabel('Time (s)')
 suptitle('First-Lick Distribution')
 
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
   
 
 % SNR plot 3 - HitRate As a function of Animal Age
%Age = 
 
%figure
%title('Hit Rate by Age')



 
 
 