old_DFF = squeeze(nanmean(OldFRA.DFF,2));
young_DFF = squeeze(nanmean(Passive.DFF,2));
species = cat(1,ones(477,1), ones(921,1)*2);

 all_DFF = cat(2,old_DFF,young_DFF);

[~,~,~,~,explained]= pca(all_DFF');
figure
plot(explained);

%%find clusters from all data 
y = tsne(all_DFF','NumPCAComponents',6);
gscatter(y(:,1),y(:,2),species)

%% find number of classes young only
figure
y = tsne(young_DFF','NumPCAComponents',2,'perplexity',30);
gscatter(y(:,1),y(:,2),Passive.Classes)