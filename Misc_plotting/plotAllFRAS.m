function  plotAllFRAS(Data)

uFreqs = unique(Data.FreqLevelOrder{:,1}); 
uLevels = unique(Data.FreqLevelOrder{:,2}) ; 


figure
for ii = 1:size(Data.DFF,3)
subplot(10,10,mod(ii-1,100)+1)
myFRA(uFreqs,uLevels', Data.df_by_level(:,ii),ii,false)
title(sprintf('neuron %d',ii))

set(gca,'CLim',[0 20])
if mod(ii,100) == 0
pause
end
end

figure
for ii = 1:size(Data.DFF,3)
subplot(10,10,mod(ii-1,100)+1)
myFRA(uFreqs,uLevels', Data.df_by_level_offset(:,ii)-Data.df_by_level(:,ii)...
    ,ii,false)
title(sprintf('neuron %d',ii))
set(gca,'CLim',[0 20])
if mod(ii,100) == 0
pause
end
end