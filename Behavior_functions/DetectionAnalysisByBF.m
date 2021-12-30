function DetectionAnalysisByBF(TN)

N_expts = size(TN,1);

BF_Neurons =[];
Neuron_norm_acc = [];

print_size = 0;

uneven_expts = [];
 
for expt = 1:N_expts
      fprintf(repmat('\b',1,print_size)) 
      print_size = fprintf('Detect: Experiment %d/%d \n',expt,N_expts);

    
    [BF_neurons_expt, BF_field] = BestFrequencyAnalysis(TN{expt,1});
    [dstats,Neuron_accuracy_expt] = DetectionAnalysisExpt(TN{expt,2});
    
    if length(BF_neurons_expt{1}) ~= length(Neuron_accuracy_expt)
         continue;
    end 
    % normalize single neuron accuracy to mean accuracy 
    Neuron_norm_acc_expt = Neuron_accuracy_expt -dstats.mean_single_neuron_prediction;
    
   if isempty(BF_Neurons) || isempty(Neuron_norm_acc)
        BF_Neurons = BF_neurons_expt{1};
        Neuron_norm_acc = Neuron_norm_acc_expt;
   else
        BF_Neurons = cat(1,BF_Neurons,BF_neurons_expt{1});
        Neuron_norm_acc = cat(1,Neuron_norm_acc,Neuron_norm_acc_expt);
   end 
     
end 
 
figure;boxplot(Neuron_norm_acc,BF_Neurons)
xlabel('Best Frequency')
ylabel('Normalized Accuracy') 
title('Single Neuron Accuracy by BF')
 
end 



