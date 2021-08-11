function Out = CompareBayes(g1,g2,GroupNames,Levels,SaveName)

% unpack
g1_time = squeeze(g1.TimeLossTotal(:,end,:));
g2_time = squeeze(g2.TimeLossTotal(:,end,:));

g1_lvl = squeeze(max(g1.TimeLossLvl(:,:,end,:)));
g2_lvl = squeeze(max(g2.TimeLossLvl(:,:,end,:)));

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

title_str = sprintf('%s%s-%s-accuracyTotal',GroupNames{1},GroupNames{2},SaveName);
title( title_str)

saveas(gcf,sprintf('%s-TonesAccuracyTiming.pdf', title_str))

[Out.all_stats,Out.main_effects]=Compare2Anova(g1_lvl,g2_lvl,GroupNames,Levels,SaveName)




