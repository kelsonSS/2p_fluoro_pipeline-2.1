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
BayesModels.Baseline_AllBehavior = [];

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
uLevels = [70,60,50];
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
    

% Create Behavioral outcomes 
Hits = Active.handles{1}.Hits(1:trials)';
Early  = Active.handles{1}.Early';
AllBehavior = Active.handles{1}.Hits +(Active.handles{1}.Early * 2) ;
AllBehavior = AllBehavior(1:trials)';


 
 for run = 1: ceil(num_reps/3)  
    
         parfor neurons = 2:num_neurons

             
             %% subscripting
             
             rand_idx = randi( size(AllData,3),neurons,1);
             mdlData = squeeze(nanmean(AllData(tone_on:tone_off,1:trials,rand_idx)));
             %% Modeling Hits
             
            
             Prediction = RunModel(mdlData,Hits);
             NumbersLossTotal(neurons) = getPredictionAccuracy(Prediction,Hits);
             %byLevel
             NumbersLossLevel(:,:,neurons) = getPredictionAccuracyByLevel(...
                                            Prediction,Hits,Levels,uLevels);
            
             %% Modeling Early                            
             Prediction = RunModel(mdlData,Early);
             NumbersLossTotal_Early(neurons) = getPredictionAccuracy(Prediction,Early);
             %byLevel
             NumbersLossLevel_Early(:,:,neurons) = getPredictionAccuracyByLevel(...
                                            Prediction,Early,Levels,uLevels)  ;      
             
             
          
             
            
             %% Modeling All behavior
             %total
             AllBehaviorPrediction =  RunModel(mdlData,AllBehavior);
             NumbersLossTotal_AllBehavior(neurons) = getPredictionAccuracy(Prediction,AllBehavior);
            
            NumbersLossLevel_AllBehavior(:,:,neurons) = getPredictionAccuracyByLevel(AllBehaviorPrediction,AllBehavior,Levels,uLevels);
            
            
            
         end
    
          
        %% Packaging  
    BayesModels.NumbersLossTotal_AllBehavior(:,expt,run) = NumbersLossTotal_AllBehavior;
    BayesModels.NumbersLossTotal_Hits(:,expt,run) = NumbersLossTotal; 
    BayesModels.NumbersLossTotal_Early(:,expt,run) = NumbersLossTotal_Early; 

    BayesModels.NumbersLossLevel_AllBehavior(:,:,expt,run) = NumbersLossLevel_AllBehavior;
    BayesModels.NumbersLossLevel_Hits(:,:,expt,run) = NumbersLossLevel;
    BayesModels.NumbersLossLevel_Early(:,:,expt,run) = NumbersLossLevel_Early;
 end
     
 
% Modeling Accuracy Over Time 
for run = 1:num_reps
    fprintf(repmat('\b',1,print_size))
    print_size = fprintf('Bayes: Expt: %d/%d Run: %d/%d - %d min  \n',expt,N_expts,run,num_reps,floor(toc/60));

    %% create Crossvalidation partitions 
   cv_Hits = cvpartition(Hits,'KFold',5);
   cv_AllBehavior = cvpartition(AllBehavior,'KFold',5);
 parfor Time = (padsize+1:90+padsize)
            
            
          
            %% subscripting
            
          %subsetting
            mdlData = AllData(Time-padsize:Time+padsize,:,:);       
            rand_idx = randi(size(mdlData,3),80,1);
           
             
            mdlData = squeeze(nanmean(mdlData(:,1:trials,rand_idx)));
            mdlData(isnan(mdlData)) = 0;
            
        
             %% Modeling Hits
             
            
             Prediction = RunModel(mdlData,Hits);
             TimeLossTotal(Time-padsize) = getPredictionAccuracy(Prediction,Hits);
             %byLevel
             TimeLossLevel(:,:,Time-padsize) = getPredictionAccuracyByLevel(...
                                            Prediction,Hits,Levels,uLevels);
              
   %% Modeling Early
             
            
             Prediction = RunModel(mdlData,Early);
             TimeLossTotal_Early(Time-padsize) = getPredictionAccuracy(Prediction,Early);
             %byLevel
             TimeLossLevel_Early(:,:,Time-padsize) = getPredictionAccuracyByLevel(...
                                            Prediction,Early,Levels,uLevels);
            
             %% Modeling All
             
            
             Prediction = RunModel(mdlData,AllBehavior)
             TimeLossTotal_AllBehavior(Time-padsize) = getPredictionAccuracy(Prediction,AllBehavior)
             %byLevel
             TimeLossLevel_AllBehavior(:,:,Time-padsize) = getPredictionAccuracyByLevel(...
                                            Prediction,AllBehavior,Levels,uLevels)
            
                      
         
 
           
        

            
 end


 BayesModels.TimeLossTotal_Hits(:,expt,run) = TimeLossTotal;
 BayesModels.TimeLossTotal_AllBehavior(:,expt,run) = TimeLossTotal_AllBehavior;
  BayesModels.TimeLossTotal_Early(:,expt,run) = TimeLossTotal_Early;
 
 BayesModels.TimeLossLevel_Hits(:,:,expt,run) = TimeLossLevel;
 BayesModels.TimeLossLevel_AllBehavior(:,:,expt,run) = TimeLossLevel_AllBehavior;
 BayesModels.TimeLossLevel_Early(:,:,expt,run) = TimeLossLevel_Early;
end

% what is model accuraccy if we guessed every trial was a hit?
BayesModels.Baseline_Hits(end+1) = sum( Hits == 1) / length(Hits);
BayesModels.Baseline_AllBehavior(end+1) = sum(AllBehavior == 1) / length(AllBehavior);

BayesModels.UsedExperimentIndex(end+1) = expt;
end 


 

end


function Prediction = RunModel(mdlData,Labels) 

 mdlData(isnan(mdlData)) = 0;
             cv = cvpartition(Labels,'KFold',5);
             mdl =  fitcnb(mdlData,Labels,'CrossVal','on','CVPartition',cv) ;
             %% Predicting
             %total
             Prediction =   mdl.kfoldPredict;


end

function acc = getPredictionAccuracy(Prediction,Labels)
    acc = sum( Prediction == Labels )...
            /  length(Labels);
%              
end 

function out = getPredictionAccuracyByLevel(Prediction,actual,Levels,uLevels)


 for lvl = 1: length(uLevels)
                 lvl_idx = Levels == uLevels(lvl);
                 if any(lvl_idx)
                  out(lvl,:) = sum( Prediction(lvl_idx) == ...
                      actual(lvl_idx) ) / length(Levels(lvl_idx));
                 else
                     out(lvl,:) = nan;
               
                 end
 
 end

end 

