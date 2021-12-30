function out = DetectionAnalysis(Active,old_idx)

if ~iscell(Active)
    Active = {Active};
end

if exist('old_idx','var');Active = Active(old_idx);end 

print_size = 0;
n_expts = length(Active);
for expt= 1:n_expts
    
   print_size = fprintf('Detect: Experiment %d/%d \n',expt,n_expts);

 temp = DetectionAnalysisExpt(Active{expt});
if expt == 1
 out = temp;
else 

out = ConcatenateStructs(out,temp);
end 


end
% are there differences between percentage of active neurons on hit v miss
% trials 
figure 
histogram(out.prc_active_hits - out.prc_active_miss,30)
title('Change in Active Neurons Hit-Miss (%)')

 


figure

% how close is animal performance to mean single neuron performance?
detect_table = struct2table(structfun(@(x) x',out,'UniformOutput',0))
figure
mdl = fitlm(detect_table , 'mean_single_neuron_prediction~animal_hit_rate')
plot(mdl)
title('all')
end









%todo- add by level 







 











