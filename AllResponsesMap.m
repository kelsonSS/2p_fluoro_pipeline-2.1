function AllResponsesMap (Raw) 

% average response over trials
DFF = Raw.DFF;

thresh = 1 ; % p <.05
Active = Raw.active.Activity;
Active = Active > thresh ;

DFF = DFF(:,:,Active); 

DFF2 = squeeze(mean(DFF,2));
Active = Raw.active.Activity;



% sort by average  max response 

[~, DFF_order] = sort(max(DFF2));
DFF2 = DFF2(:,DFF_order);

% plot average responses 
imagesc(DFF2')

[DFF2_max,DFF2_idx] = max(DFF2);

[~,DFF3_order] = sort(DFF2_idx); 

DFF3 = DFF2./ repmat(max(DFF2),[150,1]);
 
DFF3 = DFF3(:,DFF3_order); 

figure

imagesc(DFF3')

set(gca,'CLim', [0 1])

end 
 