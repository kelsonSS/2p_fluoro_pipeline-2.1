function Out = CompareBayes(g1,g2,GroupNames,Levels,SaveName)

% unpack
g1_time = squeeze(g1.TimeLossTotal(:,end,:));
g2_time = squeeze(g2.TimeLossTotal(:,end,:));

g1_neurons = squeeze(g1.ClassesTonesLossTotal(:,end,:));
g2_neurons = squeeze(g2.ClassesTonesLossTotal(:,end,:));

g1_lvls = squeeze(max(g1.TimeLossLvl(:,:,end,:)));
g2_lvls = squeeze(max(g2.TimeLossLvl(:,:,end,:)));

figure
hold on 
shadedErrorBar([],nanmean(g1_time'),nanstd(g1_time'),'b')
shadedErrorBar([],nanmean(g2_time'),nanstd(g2_time'))
h = legend(GroupNames{1},'','','',GroupNames{2})
h.String = h.String(1:4:end)
xticks([0:30:150])
xticklabels([0:1:5])
xlabel('Seconds')
ylabel('Fraction Correct')

title_str = sprintf('%s%s-%s-TonesAccuracyTiming',GroupNames{1},GroupNames{2},SaveName);
title( title_str)

saveas(gcf,sprintf('%s.pdf', title_str))


figure
hold on
shadedErrorBar([],nanmean(g1_neurons'),nanstd(g1_neurons'),'b')
shadedErrorBar([],nanmean(g2_neurons'),nanstd(g2_neurons'))
h = legend(GroupNames{1},'','','',GroupNames{2})
h.String = h.String(1:4:end)
xlabel('Neurons')
ylabel('Fraction Correct')
title(title_str)
title_str = sprintf('%s%s-%s-TonesAccuracyNeurons',GroupNames{1},GroupNames{2},SaveName);
saveas(gcf,sprintf('%s.pdf', title_str))


n_neurons = arrayfun(@num2str,[1:100],'UniformOutput',0);
n_times = arrayfun(@num2str,[1:size(g1_time,1)],'UniformOutput',0);

[Out.Neurons.all_stats,Out.Neurons.main_effects] = Compare2Anova(g1_neurons,g2_neurons,GroupNames,n_neurons,SaveName)
[Out.Time.all_stats,Out.Time.main_effects] = Compare2Anova(g1_time,g2_time,GroupNames,n_times,SaveName)
[Out.Lvls.all_stats,Out.Lvls.main_effects]=Compare2Anova(g1_lvls,g2_lvls,GroupNames,Levels,SaveName)




