function Out = PlotBayesActive(BayesModels,savename)
% creates the plots for the bayes active paper
% this function takes a BayesModel struct from the Output of
% BayesClassifierPassive and plots the resulting models 

if ~exist('savename','var')
    savename = [];
end 

if isfield(BayesModels,'BayesModels')
    BayesModels = BayesModels.BayesModels;
end 
    
%if ~isfield(BayesModels,'NumbersLossTotal')
%    error('Object must be BayesModelsObject or Contain a field named BayesModels')
%end

fn = fieldnames(BayesModels);


for model = 1:length(fn)
    
    if ~contains(fn{model},'Total')
        continue
   
    end 
        
   
    
    data = BayesModels.(fn{model});
    
    data = permute(data,[1,3,2]); % time x repeats x experiment
    data = data(4:end,:,:);
    
    
    
    
    % remove expts with no predictions 
    no_predict_idx  = squeeze(nanmean(nanmean(data))) == 0 ;
    data(:,:,no_predict_idx  ) =[];
    % smooth data
    
    data = smoothdata(data,1,'movmean',5);
    
%     data = BayesModels.(fn{model});
%     data = permute(data(4:end,:,:),[1,3,2])
%     
     B_Vec = repmat(nanmean(data(1:30,:,:)),[size(data,1),1,1]);
     %Vec_DFF = (data -B_Vec)./B_Vec * 100;
     DF = data - B_Vec;
     %Vec_DFF_mu = squeeze(nanmean(Vec_DFF,2));
     Vec_DF_mu = squeeze(nanmean(DF,2));
     
     
     %% begin plotting 
     fig_name = ['Bayes-' fn{model} ];
     figure
     set(gcf,'Position', [1 1 800 1200])
     hold on
     suptitle(strrep(fn{model},'_','-'))
     
     
     subplot(2,2,1)
     hold on
     plot(squeeze(nanmean(data,2)))
     plot(squeeze(nanmean(nanmean(data,2),3)),'k','LineWidth',3)
     xticks(0:30:91)
     xticklabels(0:1:4)
     ylim([0 1])
     xlabel('time (s)')
     ylabel('delta accuracy')
     
     
     % average delta accuracy
     subplot(2,2,2)
     BayesFig(Vec_DF_mu)
     
     % Seperate by max  delta prediction
     max_prediction = squeeze(max(nanmean(data,2)));
     min_prediction = squeeze(min(nanmean(data,2)));
     max_delta_prediction = max_prediction - min_prediction;
     neural_decode_idx = max_delta_prediction  >= .1 ;
     
     n_good_decode = sum(neural_decode_idx);
     n_bad_decode = sum(~neural_decode_idx);
     total = length(neural_decode_idx);
     % plotting good
     subplot(2,2,3)
     title(sprintf('Neural Decoding: \n %d of %d',n_good_decode,total))
     if n_good_decode
     BayesFig(DF(:,:,neural_decode_idx));
     end 
    title(sprintf('Neural Decoding: \n %d of %d',n_good_decode,total))
     % Plotting bad 
     subplot(2,2,4)
     if n_bad_decode
     BayesFig(DF(:,:,~neural_decode_idx));
     end 
      title(sprintf('No Decoding:\n %d of %d',n_bad_decode,total))
  
    %% saving figure 
    if savename
        outpath = [savename '-' fig_name, '.pdf'];
        saveas(gcf,outpath) 
        close(gcf)
    end 
    
     
     
     %% packaging output
     Out.(fn{model}).MaxPrediction = max_prediction;
     Out.(fn{model}).MaxDeltaPrediction = max_delta_prediction(:);
     
     
        
              
end




function BayesFig(data)
% plots and styles the figures

data = squeeze(data);

if ndims(data) == 3
    data = reshape(data,size(data,1),[]);
end

% data should now be 2d steps x reps matrix

shadedErrorBar([],nanmean(data,2),nanstd(data,[],2))
ylim([-.1,.3]) 
xticks(0:30:91)
xticklabels(0:1:4)
xlabel('time (s)')
ylabel('delta accuracy')





  





