function PlotFirstResponse(TNBehavior,SaveName)


first_licks = cell2mat(...
    cellfun(@(x)  x.handles{1}.FirstResponse(:,1),...
             TNBehavior,'UniformOutput',0));

figure
histogram(first_licks,0:.1:5,'Normalization','probability')
hold on

ylim([0 .2])
ax = gca;
plot([1.2 1.2], [ax.YLim(1) ax.YLim(2)],'k--')
plot([2.5 2.5], [ax.YLim(1) ax.YLim(2)],'k--')
ticks = yticks;
yticklabels(ticks * 100);
ylabel('%')
xlabel('Time (S)')




if exist('SaveName','var')
    title(SaveName,'Interpreter','none')
    saveas(gcf,sprintf('%s.pdf',SaveName))
end 