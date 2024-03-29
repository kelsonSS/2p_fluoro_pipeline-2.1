function PlotBayes(BayesModels,savepath)

% this function takes a BayesModel struct from the Output of
% BayesClassifierPassive and plots the resulting models 

if ~exist('savepath','var')
    savepath = []
end 

if isfield(BayesModels,'BayesModels')
    BayesModels = BayesModels.BayesModels;
end 
    
%if ~isfield(BayesModels,'NumbersLossTotal')
%    error('Object must be BayesModelsObject or Contain a field named BayesModels')
%end

fn = fieldnames(BayesModels);


for model = 1:length(fn)
    

        
    
    data = BayesModels.(fn{model});
    if ndims(data) == 2
        continue
    end
    if ndims(data) == 3 
    data = reshape(data,size(data,1),1,size(data,2),size(data,3) );    
    end 
    
%     data = BayesModels.(fn{model});
%     data = permute(data(4:end,:,:),[1,3,2])
%     
%     B_Vec = repmat(nanmean(data(1:30,:,:)),[90,1,1]);
%     Vec_DFF = (data -B_Vec)./B_Vec * 100;
%     DF = data - B_Vec;
%     Vec_DFF_mu = squeeze(nanmean(Vec_DFF,2))
%     Vec_DF_mu = squeeze(nanmean(DF,2))
%     plotAllDFFs(data)
    % usually 1 but could be larger if more expts run or if done on per
    % animal basis
    for expt = 1:size(data,3)
        
        figure
         hold on 
         fig_name = [fn{model} ' Experiment: ' num2str(expt) ];
         
        
         for lvl = 1:size(data,2)
            BayesFig(data(:,lvl,expt,:))
         end 
    
     saveStyleBayesFig(fig_name,savepath)         
              
    end 
    
    
    
end 




function BayesFig(data)
% plots and styles the figures

data = squeeze(data);

% data should now be 2d steps x reps matrix

shadedErrorBar([],nanmean(data,2),nanstd(data,[],2))


function saveStyleBayesFig(fig_name,savepath)


x_lim = xlim;
xlim([2 x_lim(2)] )
ylabel('Fraction Correct')
title(fig_name)
if ~isempty(savepath) 
    fig_name =  matlab.lang.makeValidName(fig_name)
   %outpath = [savepath '.ps'];

     outpath = fullfile(savepath,[fig_name, '.pdf']);
    print(outpath,'-dpdf','-bestfit')
   % print(outpath, '-dpsc', '-append') 
   pause(1)
    close(gcf)
end






