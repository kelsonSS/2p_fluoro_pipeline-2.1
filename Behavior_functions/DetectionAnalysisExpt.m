function [out,single_neuron_accuracy,N_activity] = DetectionAnalysisExpt(A)
% given an active dataset (A) created with fluoro_to_table
% returns 
% out - various statistics of the behavioral performance of the animal 
% single_neuron_predictions - a matrix of which neurons were active on each
%                             trial 

out = [];
n_neurons = size(A.DFF,3);
n_trials = size(A.DFF,2);
handles = A.handles{1};

if length(handles.Hits) < n_trials
    n_trials = length(handles.Hits)
end 

Hits = logical(handles.Hits(1:n_trials));
Early = logical(handles.Early(1:n_trials)); 
Miss = logical(handles.Miss(1:n_trials));



for trial_idx = 1:n_trials
for neuron_idx = 1:n_neurons

N_activity(neuron_idx,trial_idx) = ttest2(A.DFF(1:30,trial_idx,neuron_idx),...
                                 A.DFF(31:60,trial_idx,neuron_idx) );
 end 
end 



single_neuron_accuracy = sum(N_activity == Hits,2) ./ n_trials;

out.single_neuron_accuracy = {single_neuron_accuracy};
out.mean_single_neuron_prediction = mean(single_neuron_accuracy);
out.max_single_neuron_prediction = max(single_neuron_accuracy);
out.animal_hit_rate = sum(Hits)/length(Hits);

total_active = sum(N_activity);
out.prc_active = mean(total_active) / n_neurons;
out.prc_active_hits = mean(total_active(Hits)) / n_neurons ;
out.prc_active_miss = mean(total_active(Miss)) / n_neurons;
out.prc_active_early = mean(total_active(Early)) / n_neurons;


