function SingleNeuronModelLevels(TNFullBehaviorTable,SaveName)

if ~exist('SaveName','var')
    SaveName = []
end 

levels = [20,10,0]
h =figure
set(h,'Position',[1280 300 400 700])
n_levels = length(levels)
for lvl = 1:n_levels

modelspec = sprintf('mean_single_neuron_prediction~dPrimeLevel_%d',levels(lvl));
Model = fitlm(TNFullBehaviorTable,modelspec);
subplot(n_levels,1,lvl)
plot(Model)
xlim([0 3])
ylim([.4 .8])
end 
if SaveName
    saveas(gcf, [SaveName '.pdf'])
end


