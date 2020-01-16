% average response over trials
DFF  = out2.DFF;
DFF2 = squeeze(mean(out2.DFF,2));

% sort by average  max response 

[~, DFF_order] = sort(max(DFF2));
DFF2 = DFF2(:,DFF_order);

% plot average responses 
imagesc(DFF2')

[DFF2_max,DFF2_idx] = max(DFF2);

[~,DFF3_order] = sort(DFF2_idx); 

DFF3 = DFF2./ repmat(max(DFF2),[150,1]);
 
DFF3 = DFF3(:,DFF3_order); 