function FramesToSeconds 
ax = gcf;
figure(ax) 
xticks([30,60,90,])
xticklabels({'1','2','3'})
xlabel('Time (sec)')
ylabel('Neurons')
caxis([0 1])
end 

