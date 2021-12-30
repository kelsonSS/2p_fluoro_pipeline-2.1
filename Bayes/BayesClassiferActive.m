function BayesModels = BayesClassiferActive(Active_All)

% this function takes in the ActiveData and creates a number of
% classifiers and looks at their classification performance based on number
% of neurons, time during the trial , and Neuron Class



%if ~exist('Savepath','var')
%Savepath = fullfile('Z:\Kelson', Savename);
%end 
N_expts = length(Active_All);
BayesModels.UsedExperimentIndex =[];
BayesModels.Baseline_Hits = [];
BayesModels.Baseline_HitsEarly = [];

padsize = 3;
print_size = 0;
tic;
for expt = 1:N_expts
Active = Active_All{expt};

tone_on = 31;
tone_off = 60;


num_neurons = 80;

num_reps = 20;

Levels = Active.FreqLevelOrder{:,2};
trials = length(Levels);
uLevels = sort(unique(Levels),'descend');
clean_idx = Active.Clean_idx;

if trials < 20 
    continue
end 

if isfield(Active,'Combined_Classes')
    all_classes = Active.Combined_Classes;
    num_classes = max(all_classes);
elseif ~(isfield(Active,'Classes'))
    all_classes = ones(size(Active.DFF,3),1);
elseif iscell(Active.Classes)
    all_classes = Active.Class_idx;
    classes =  Active.Classes;  
else 
    all_classes = Active.Classes;
end 

%% Choose Responsive Neuron Critieria (pick one)

% neurons that significantly respond to at least one stimulus( less restrictive)

%active_idx = Active.responsive;

% Active neurons were responsive to sound overall ( more restrictive)

active_idx = Active.active{:,2}>0;

if sum(active_idx & clean_idx) < num_neurons
    active_idx = Active.responsive;
end 

in_idx =  active_idx & clean_idx ;
          
        AllData = Active.DFF_Z(:,:,in_idx);          
        if size(AllData,3) < num_neurons  
            continue
        end
    

% detect_Hits 
Hits = Active.handles{1}.Hits(1:trials)';
HitsEarly = Active.handles{1}.Hits +(Active.handles{1}.Early * 2) ;
HitsEarly = HitsEarly(1:trials)';


 
 for run = 1: ceil(num_reps/3)  

          
         parfor neurons = 2:num_neurons
             
             
             %% subscripting
             
             rand_idx = randi( size(AllData,3),neurons,1);
             %% Modeling
             mdlData = squeeze(nanmean(AllData(tone_on:tone_off,1:trials,rand_idx)));
             mdlData(isnan(mdlData)) = 0;
             cv = cvpartition(Hits,'KFold',5);
             mdl =  fitcnb(mdlData,Hits,'CrossVal','on','CVPartition',cv) ;
             %% Predicting
             %total
             Prediction =   mdl.kfoldPredict;
             NumbersLossTotal(neurons) = sum( Prediction == Hits )...
                 /length(Hits);
%              % by level
%              for lvl = 1: length(uLevels)
%                  lvl_idx = Levels == uLevels(lvl);
%                  Active.BayesModels.NumbersLossLvl(lvl,neurons,expt,run) = sum( Prediction(lvl_idx) == ...
%                      Freqs(lvl_idx) )/length(Freqs(lvl_idx));
%              end
%               mdlData(isnan(mdlData)) = 0;
             cv = cvpartition(HitsEarly,'KFold',5);
             mdl =  fitcnb(mdlData,HitsEarly,'CrossVal','on','CVPartition',cv) ;
             %% Predicting
             %total
             Prediction =  mdl.kfoldPredict;
             NumbersLossTotal_HitsEarly(neurons) = sum( Prediction == HitsEarly )...
                 /length(HitsEarly);
            
             
         end
    BayesModels.NumbersLossTotal_Hits(:,expt,run) = NumbersLossTotal;
    BayesModels.NumbersLossTotal_HitsEarly(:,expt,run) = NumbersLossTotal_HitsEarly; 
 end
     
%  
% 
% for run = 1:num_reps
%     fprintf(repmat('\b',1,print_size))
%     print_size = fprintf('Bayes: Expt: %d/%d Run: %d/%d - %d min  \n',expt,N_expts,run,num_reps,floor(toc/60));
% 
%     %% create Crossvalidation partitions 
%    cv_Hits = cvpartition(Hits,'KFold',5);
%    cv_HitsEarly = cvpartition(HitsEarly,'KFold',5);
%  parfor Time = (padsize+1:90+padsize)
%             
%             
%           
%             %% subscripting
%             
%           %subsetting
%             mdlData = AllData(Time-padsize:Time+padsize,:,:);       
%             rand_idx = randi(size(mdlData,3),80,1);
%             
%             
%             %% Modeling Hits
%              mdlData = squeeze(nanmean(mdlData(:,1:trials,rand_idx)));
%              mdlData(isnan(mdlData)) = 0;
%             
%             
%              
%              mdl =  fitcnb(mdlData,Hits,'CrossVal','on','CVPartition',cv_Hits) ;
%  
%             %% Predicting
%             
%             %total
%             Prediction =  mdl.kfoldPredict;
%             %BayesModels.TimeLossTotal(Time,expt,run) = sum( Prediction == Hits )...
%             run_TimeLossTotal(Time-padsize) = sum( Prediction == Hits )...
%                 /length(Hits);
%           
%             %% Modeling Hits & Early Trials
%              
%              mdl =  fitcnb(mdlData,HitsEarly,'CrossVal','on','CVPartition',cv_HitsEarly) ;
%  
%             %% Predicting
%             
%             %total
%             Prediction =  mdl.kfoldPredict;
%             %BayesModels.TimeLossTotal_HitsEarly(Time,expt,run) = sum( Prediction == HitsEarly )...
%             run_TimeLossTotal_HitsEarly(Time-padsize) = sum( Prediction == HitsEarly )...
%             /length(HitsEarly);
% 
%             
%  end
%   BayesModels.TimeLossTotal_Hits(:,expt,run) = run_TimeLossTotal;
%  BayesModels.TimeLossTotal_HitsEarly(:,expt,run) = run_TimeLossTotal_HitsEarly;
% end

% what is model accuraccy if we guessed every trial was a hit?
BayesModels.Baseline_Hits(end+1) = sum( Hits == 1) / length(Hits);
BayesModels.Baseline_HitsEarly(end+1) = sum(HitsEarly == 1) / length(HitsEarly);

BayesModels.UsedExperimentIndex(end+1) = expt;
end 


 

end


