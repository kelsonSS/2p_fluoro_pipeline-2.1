function DF = SubsetDFF(DF,idx)
     DF.DFF = DF.DFF(:,idx,:);
     DF.DFF_Z = DF.DFF_Z(:,idx,:);
     DF.DFF_norm = DF.DFF_norm(:,idx,:);
     DF.FreqLevelOrder = DF.FreqLevelOrder(idx,:);


end