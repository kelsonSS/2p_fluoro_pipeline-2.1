function accuracy = SingleNeuronPredictions(TNBehavior)

accuracy = [];

n_expts = length(TNBehavior)
for expt_idx = 1:n_expts
    fprintf('expt: %d of %d\n', expt_idx,n_expts)
    expt = TNBehavior{expt_idx};
    
    hits = expt.handles{1}.Hits;
    
    neuron_accuracy_expt = getAccuracyOverTime(expt,hits);
    
    accuracy = cat(2,accuracy, neuron_accuracy_expt);
    
end


function acc = getAccuracyOverTime(expt,hits)

DF = expt.DFF(:,:,logical(expt.responsive));

activity = NeuronActivityOverTime(DF);

acc = activity == hits;

acc = squeeze(nanmean(acc,2));





function active = NeuronActivityOverTime(DF)
    window_size = 5;
    active = zeros(size(DF));
    n_frames = size(DF,1);
    n_trials = size(DF,2);
    n_neurons = size(DF,3);
    padding = zeros(window_size,n_trials,n_neurons);
    DF = cat(1,padding,DF);
    DF = cat(1,DF,padding);
    for time = window_size:n_frames+window_size
        for trial = 1:n_trials
            for neuron = 1:n_neurons
                baseline = DF(window_size:2*window_size,trial,neuron);
                active(time,trial,neuron) = ttest2(baseline,...
                                   DF(time:time+window_size,trial,neuron));
            end 
        end 
    end 
        
   
    



    