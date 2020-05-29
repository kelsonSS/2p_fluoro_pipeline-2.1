function  plotAllFRAS(Data,sig)

uFreqs = unique(Data.FreqLevelOrder{:,1}); 
uLevels = unique(Data.FreqLevelOrder{:,2}) ; 

active_idx = find(Data.active{:,2} > 1 );


figure
for ii = 1:length(active_idx)
subplot(10,10,mod(ii-1,100)+1)
Neuron = active_idx(ii);
in = Data.df_by_level(:,:,Neuron);
end
if sig

myFRA(uFreqs,uLevels', in ,Neuron,false)
title(sprintf('neuron %d',ii))

%set(gca,'CLim',[0 20])
 if mod(ii,100) == 0
 pause
 end
end
