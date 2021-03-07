function  plotAllFRAS(Data,sig)

if ~exist('sig','var')
    sig = 0
end 
uFreqs = unique(Data.FreqLevelOrder{:,1}); 
uLevels = unique(Data.FreqLevelOrder{:,2}) ; 

active_idx = find(Data.active{:,2} > sig );

cell_lengths = cellfun(@(x) size(x,3),Data.df_by_level);
figure
for ii = 1:length(active_idx)
subplot(10,10,mod(ii-1,100)+1)
Neuron = active_idx(ii);
[cell_num,nn] = findLocation(Neuron,cell_lengths);
in = Data.df_by_level{cell_num}(:,:,nn);


myFRA(uFreqs,uLevels', in ,Neuron,false)
title(sprintf('neuron %d',ii))

%set(gca,'CLim',[0 20])
 if mod(ii,100) == 0
 pause
 end
end

function [cell_num,Neuron] = findLocation(Neuron,Cell_lengths) 

cell_num = 1; 
while Neuron > Cell_lengths(cell_num)
    Neuron = Neuron - Cell_lengths(cell_num);
    cell_num = cell_num +1;
end 

