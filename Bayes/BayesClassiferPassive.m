function  Loss = BayesClassiferPassive(Passive,Savename,Savepath)

% this function takes in the passiveData and creates a number of
% classifiers and looks at their classification performance based on number
% of neurons, time during the trial , and Neuron Class

% Savename - the File you want to save to. should be the same name as the
% DF 
% init


if ~exist('Savepath','var')
Savepath = fullfile('\\vault3\Data\Kelson', Savename);
end 

handles = Passive.handles{1};
tone_on = handles.PreStimSilence * 30 ;
tone_off = handles.PrimaryDuration * 30 + tone_on; 


if length(handles.BackgroundNoise) == 3 && handles.BackgroundNoise(1) ~=99
    noise_on = handles.BackgroundNoise(2) * 30;
    noise_off = handles.BackgroundNoise(3)* 30;
end 

num_classes = 4;
num_neurons = 100;
num_types = 2;
num_reps = 10;
classes = {'Tones','Noise','Offset'};

Freqs = Passive.FreqLevelOrder{:,1};
Levels = Passive.FreqLevelOrder{:,2};
trials = length(Levels);
uLevels = sort(unique(Levels),'descend');
bad_idx = ~Passive.Clean_idx;
if ~(isfield(Passive,'Classes'))
    all_classes = ones(size(Passive.DFF,3),1);
elseif iscell(Passive.Classes)
    all_classes = Passive.Class_idx;
    classes =  Passive.Classes;
    
else 
    all_classes = Passive.Classes;
end 

try
    expt_list = Passive.experiment_list;
    active_idx = Passive.active{:,2} >0 ;
    
catch
    expt_list = Passive.Experiment_list % for compatibility with old files
     active_idx = Passive.Active{:,2} >0 ;
end 



Levels2 = Levels ;
Levels2(Levels2 ~= 99) = 1;
Levels2(Levels2 == 99) = 0 ;
uLevels2 = unique(Levels2);



% test the number of neurons that are needed for classifications at different SNRS
tic;
for expt = 1:size(Passive.DataDirs)
    disp([' Numbers: expt #' num2str(expt) '/' num2str(size(Passive.DataDirs))  ': ' ...
        num2str(floor( toc / 60  )) ' min elapsed' ])
    for run = 1:ceil( num_reps / length(Passive.DataDirs) )
        %% indexing
        expt_idx = expt_list == expt;
       
        in_idx = expt_idx & active_idx & ~ bad_idx ;
        
        AllData = Passive.DFF_Z(:,:,in_idx);
        
        if size(AllData,3) < 20
            continue
            end
        
        for neurons = 2:num_neurons
            
            
            %% subscripting
            
            rand_idx = randi( size(AllData,3),neurons,1);
            %% Modeling
            mdlData = squeeze(nanmean(AllData(tone_on:tone_off,1:trials,rand_idx)));
            mdlData(isnan(mdlData)) = 0;
            mdl =  fitcnb(mdlData,Freqs,'CrossVal','on') ;
            %% Predicting
            %total
            Prediction =   mdl.kfoldPredict;
            Passive.BayesModels.NumbersLossTotal(neurons,expt,run) = sum( Prediction == Freqs )...
                /length(Freqs);
            % by level
            for lvl = 1: length(uLevels)
                lvl_idx = Levels == uLevels(lvl);
                Passive.BayesModels.NumbersLossLvl(lvl,neurons,expt,run) = sum( Prediction(lvl_idx) == ...
                    Freqs(lvl_idx) )/length(Freqs(lvl_idx));
            end
            
            
            
        end
    end
end
Passive.BayesModels.NumbersLossLvl = permute(Passive.BayesModels.NumbersLossLvl,[2,1,3,4]);

save(Savepath,'Passive','-v7.3') 

%Test how classes of neurons encode Tone  information
for run = 1:num_reps
    
     disp([' Classes Tones: Run #' num2str(run) '/' num2str(num_reps)  ': ' ...
        num2str(floor( toc / 60  )) ' min elapsed' ])
    
    for class = 1:length(classes) +1
        
        %% indexing
        class_idx = all_classes == class;                         
        if class == length(classes)+1       
            class_idx = ones(size(class_idx));   
        end             
        
      
        
        in_idx = class_idx & active_idx & ~ bad_idx;
        AllData = Passive.DFF_Z(tone_on:tone_off,:,in_idx);
        
         if size(AllData,3) < 20
            continue
            end
        for neurons = 2:num_neurons
            
            
            
            %% subscripting
            rand_idx = randi( size(AllData,3),neurons,1);
            %% Modeling
            mdlData = squeeze(nanmean(AllData(:,1:trials,rand_idx)));
             mdlData(isnan(mdlData)) = 0;
            mdl =  fitcnb(mdlData,Freqs,'Crossval','on') ;
                       
            %% Predicting
            %total
            Prediction = mdl.kfoldPredict;
            Passive.BayesModels.ClassesTonesLossTotal(neurons,class,run) = sum( Prediction == Freqs )...
                /length(Freqs);
            % by level
            
            for lvl = 1: length(uLevels)
                lvl_idx = Levels == uLevels(lvl);
                Passive.BayesModels.ClassesTonesLossLvl(lvl,neurons,class,run) = sum( Prediction(lvl_idx) == ...
                    Freqs(lvl_idx) )/length(Freqs(lvl_idx));
            end
            
            
        end
    end
end
Passive.BayesModels.ClassesTonesLossLvl = permute(Passive.BayesModels.ClassesTonesLossLvl,[2,1,3,4]);

% save(Savepath,'Passive','-v7.3') 



% %% Test How the same Classes encode Noise Information
%
%
for run = 1:num_reps
    
     disp([' Classes Noise: Run #' num2str(run) '/' num2str(num_reps)  ': ' ...
        num2str(floor( toc / 60  )) ' min elapsed' ])
    
    for neurons = 2:num_neurons
        for class = 1:length(classes)+1
            
            %% indexing
           
                   
            class_idx = all_classes == class;
            
             if class == length(classes)+1
                 class_idx = ones(size(class_idx));
             end 
           
            
            in_idx = class_idx & active_idx & ~ bad_idx;
            
            %% subscripting
            AllData = Passive.DFF_Z(noise_on:tone_on,:,in_idx);
             if size(AllData,3) < 20
            continue
            end
            rand_idx = randi( size(AllData,3),neurons,1);
            %% Modeling
            mdlData = squeeze(nanmean(AllData(:,1:trials,rand_idx)));
             mdlData(isnan(mdlData)) = 0;
            mdl =  fitcnb(mdlData,Levels2,'Crossval','on') ;
            %% Predicting
            %total
            Prediction =  mdl.kfoldPredict;
            Passive.BayesModels.ClassesNoiseLossTotal(neurons,class,run) = sum( Prediction == Levels2 )...
                /length(Levels2);
            % by level
           
            % by level- reshaping into reps x lvls x Freqs
            Prediction_by_Lvl = reshape(Prediction,10,4,[]);
            Prediction_by_Lvl = permute(Prediction_by_Lvl,[3,1,2]);
            
            Levels_by_Lvl = reshape(Levels2,10,4,[]);
            Levels_by_Lvl = permute(Levels_by_Lvl,[3,1,2]);
          
            Passive.BayesModels.ClassesNoiseLossLvl(neurons,:,class,run) =...
                squeeze(sum(sum(...
                                 Prediction_by_Lvl == Levels_by_Lvl))...
                                 / (8*10) );
            
            
        end
    end
end
 save(Savepath,'Passive','-v7.3') 


%%% Classes information over time
%% Tones Over time 

LossLvl = nan(size(Passive.DFF,1),length(uLevels),num_classes,num_reps);
tic;
frame_window = 5;
padsize = floor(frame_window/2)  ;

AllData = Passive.DFF_Z;

zData = AllData;

for run = 1:num_reps
    disp(['Time: Run #' num2str(run) '/' num2str(num_reps)  ': ' ...
        num2str(floor( toc / 60  )) ' min elapsed' ])
    for class = 1:length(classes) +1
        for Time = padsize+1:150+padsize
            
            
            %% indexing
            class_idx = all_classes == class;         
            if class == length(classes)+1                
                class_idx = ones(size(class_idx));         
            end             
            
          
            
        
            
            in_idx = class_idx & active_idx & ~ bad_idx  ;
            
            %% subscripting
            % AllData = Passive.DFF_Z;
            
            %             AllData_size = size(AllData);
            %
            %
            %              % shuffling
            %             AllData = AllData(randperm(AllData_size(1)),...
            %                   randperm(AllData_size(2)),...
            %                   randperm(AllData_size(3) ));
            %
            %            AllData(isnan(AllData)) = rand(sum(sum(sum(isnan(AllData)))),1 );
            % selection
            AllData = zData(:,:,in_idx);
            
            % padding
            padding =  nan(padsize,size(AllData,2),size(AllData,3))  ;
            AllData = [ padding;AllData;padding];
            AllData = AllData(Time-padsize:Time+padsize,:,:);
            
             if size(AllData,3) < 20
            continue
            end
            
            
            rand_idx = randi(size(AllData,3),80,1);
            
            
            %% Modeling
            mdlData = squeeze(nanmean(AllData(:,1:trials,rand_idx)));
             mdlData(isnan(mdlData)) = 0;
            
            
            mdl =  fitcnb(mdlData,Freqs,'Crossval','on') ;
            %% Predicting
            
            %total
            Prediction =  mdl.kfoldPredict;
            Passive.BayesModels.TimeLossTotal(Time,class,run) = sum( Prediction == Freqs )...
                /length(Freqs);
            PredictionsTotal (:,Time,class,run)= Prediction;
            
            % by level- reshaping into reps x lvls x Freqs
            Prediction_by_Lvl = reshape(Prediction,10,4,[]);
            Prediction_by_Lvl = permute(Prediction_by_Lvl,[3,1,2]);
            
            Freqs_by_Lvl = reshape(Freqs,10,4,[]);
            Freqs_by_Lvl = permute(Freqs_by_Lvl,[3,1,2]);
            Passive.BayesModels.TimeLossLvl(Time,:,class,run) = squeeze(sum(sum(...
                Prediction_by_Lvl == Freqs_by_Lvl))...
                / (8*10) );
            %
            %         by level- validation
            %          for lvl = 1: length(uLevels)
            %                 lvl_idx = Freqs == uLevels(lvl);
            %                 LossLvl(Time,lvl,class,run) = sum( Prediction(lvl_idx) == ...
            %                     Freqs(lvl_idx) )/length(Freqs(lvl_idx));
            %             end
            %
            
        end
    end
end

save(Savepath,'Passive','-v7.3') 


%% Noise Over time 

LossLvl = nan(size(Passive.DFF,1),length(uLevels),num_classes,num_reps);
tic;
frame_window = 5;
padsize = floor(frame_window/2)  ;

AllData = Passive.DFF_Z;

zData = AllData;

for run = 1:num_reps
    disp(['Time-Noise: Run #' num2str(run) '/' num2str(num_reps)  ': ' ...
        num2str(floor( toc / 60  )) ' min elapsed' ])
    for class = 1:length(classes) +1
        for Time = padsize+1:150+padsize
            
            
            %% indexing
            class_idx = all_classes == class;  
            if class == length(classes)+1      
                class_idx = ones(size(class_idx));      
            end             
            
           
        
            
            in_idx = class_idx & active_idx & ~ bad_idx  ;
            
            %% subscripting
            % AllData = Passive.DFF_Z;
            
            %             AllData_size = size(AllData);
            %
            %
            %              % shuffling
            %             AllData = AllData(randperm(AllData_size(1)),...
            %                   randperm(AllData_size(2)),...
            %                   randperm(AllData_size(3) ));
            %
            %            AllData(isnan(AllData)) = rand(sum(sum(sum(isnan(AllData)))),1 );
            % selection
            AllData = zData(:,:,in_idx);
            
             if size(AllData,3) < 20
            continue
            end
            
            % padding
            padding =  nan(padsize,size(AllData,2),size(AllData,3))  ;
            AllData = [ padding;AllData;padding];
            AllData = AllData(Time-padsize:Time+padsize,:,:);
            
            
            
            rand_idx = randi(size(AllData,3),80,1);
            
            
            %% Modeling
            mdlData = squeeze(nanmean(AllData(:,1:trials,rand_idx)));
             mdlData(isnan(mdlData)) = 0;
            
            
            mdl =  fitcnb(mdlData,Levels2,'Crossval','on') ;
            %% Predicting
            
            %total
            Prediction =  mdl.kfoldPredict;
            Passive.BayesModels.TimeLossNoiseTotal(Time,class,run) = sum( Prediction == Levels2 )...
                /length(Freqs);
            PredictionsTotal (:,Time,class,run)= Prediction;
            
            % by level- reshaping into reps x lvls x Freqs
            Prediction_by_Lvl = reshape(Prediction,10,4,[]);
            Prediction_by_Lvl = permute(Prediction_by_Lvl,[3,1,2]);
            
            Levels_by_Lvl = reshape(Levels2,10,4,[]);
            Levels_by_Lvl = permute(Levels_by_Lvl,[3,1,2]);
            Passive.BayesModels.TimeLossNoiseLvl(Time,:,class,run) = squeeze(sum(sum(...
                Prediction_by_Lvl == Levels_by_Lvl))...
                / (8*10) );
            %
            %         by level- validation
            %          for lvl = 1: length(uLevels)
            %                 lvl_idx = Freqs == uLevels(lvl);
            %                 LossLvl(Time,lvl,class,run) = sum( Prediction(lvl_idx) == ...
            %                     Freqs(lvl_idx) )/length(Freqs(lvl_idx));
            %             end
            %
            
        end
    end
end

save(Savepath,'Passive','-v7.3') 


% 
% 
% 
% for run = 1:num_reps
%     
%      disp(['Tone Classes Tones: Run #' num2str(run) '/' num2str(num_reps)  ': ' ...
%         num2str(floor( toc / 60  )) ' min elapsed' ])
%     
%     
%         %% indexing
%         class_idx = (all_classes == 2 | all_classes == 3 ) ;
%         
%       
%         
%         in_idx = class_idx & active_idx & ~ bad_idx;
%         AllData = Passive.DFF_Z(tone_on:tone_off,:,in_idx);
%         
%         for neurons = 2:num_neurons
%             
%             
%             
%             %% subscripting
%             rand_idx = randi( size(AllData,3),neurons,1);
%             %% Modeling
%             mdlData = squeeze(nanmean(AllData(:,1:trials,rand_idx)));
%             mdl =  fitcnb(mdlData,Freqs,'Crossval','on') ;
%             %% Predicting
%             %total
%             Prediction = mdl.kfoldPredict;
%             Passive.BayesModels.MultiClassesTonesLossTotal(neurons,class,run) = sum( Prediction == Freqs )...
%                 /length(Freqs);
%             % by level
%             
%             for lvl = 1: length(uLevels)
%                 lvl_idx = Levels == uLevels(lvl);
%                 Passive.BayesModels.ToneClassesTonesLosslvl(lvl,neurons,class,run) = sum( Prediction(lvl_idx) == ...
%                     Freqs(lvl_idx) )/length(Freqs(lvl_idx));
%             end
%             
%             
%         end
% end
% 
%   Passive.BayesModels.ToneClassesTonesLosslvl =  permute(...
%       Passive.BayesModels.ToneClassesTonesLosslvl,[ 2,1,3,4]);
% 
% save(Savepath,'Passive','-v7.3') 





% theoretical minimum
x =  randi(8,320,2,80,10000);
y = squeeze( sum(x(:,1,:,:) == x(:,2,:,:))./320);



% %% plotting
% % all
% temp = reshape(LossTotal,size(LossTotal,1),[]);
% figure
% shadedErrorBar([],mean(temp,2),std(temp,[],2)/sqrt(60 )* 1.96 )
% hold on
% shadedErrorBar([],mean(y,2),std(y,[],2)/sqrt(80)*1.96)
%
% % by Level
% temp = reshape(LossLvl,length(uLevels),size(LossLvl,2),[]);
% temp = permute(temp,[2,1,3]);
% figure
% hold on
% colors = {'b','r','g','k'}
% for lvl = 1:length(uLevels)
% shadedErrorBar([],mean(temp(:,lvl,:),3), std(temp(:,lvl,:),[],3)...
%                                          /sqrt(size(temp,3 ))* 1.96,colors(lvl)) %95-CI
% end
% shadedErrorBar([],mean(y,2),std(y,[],2)/sqrt(80)*1.96)


% Dans suggestion to add if merging classes can give more information than
% using one class alone ( AKA if P(Tone-On + Tone-OFF) > P(Tone-On)



toc;
%TODO
%
%
% % Test in vs out of network neurons are better at classification
% types = {'All','Tones','Noise','Offset'}
% for expt = 1:size(Passive.DataDirs)
%     for type = 1:3
%     idx = expt_list == expt;
%     AllData = Passive.DFF(:,:,expt);
%     in_idx = Passive.GrangerCellIDS.(types{type}){expt}
%      AllData,FreqLevelOrder(:,2))
%

% test timing of signals
