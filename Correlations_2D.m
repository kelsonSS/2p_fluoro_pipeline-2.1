function CorrsBetween =  Correlations_2D(Passive,use_clusters,savepath)

if ~exist('Passive','var')
    load('//vault3/data/kelson/Passive/')
end

if ~exist('use_clusters','var');use_clusters = 0;end
if ~exist('savepath','var');savepath = [];end

if use_clusters
    
    names = Passive.Classes;
    num_clust = length(names);
    
    for C = 1:num_clust
        indices{C} = Passive.Class_idx == C;
    end 
    
    
    
else 
   error('using handpicked clusters disabled') 
% % Previous Code for hand picked clusters - legacy code  
% names = {'Tone_on','Noise_on','Tone_off'};
% DFF_mu = squeeze(nanmean(Passive.DFF,2));
% [~,onsets]  = max(DFF_mu);
% 
% Noise_idx = onsets > 30 & onsets <= 60;
% Tone_idx  = onsets > 60 & onsets <= 90;
% Offset_idx = onsets > 90 & onsets<= 120;
% 
% 
% indices = {Tone_idx,Noise_idx,Offset_idx};
end 
% get triangular matrix as indices 

 num_classes = length(indices);
 idx = nchoosek( 1:num_classes,2 );
 
 % append main diagonal of triangular matrix
 idx = cat(1,idx, [1:num_classes;1:num_classes]');
 
 
 
    
for index_num =  1:length(idx);
   
     if exist('progress_text_length','var')
       fprintf(repmat('\b',1,progress_text_length)) % clear previous line
   end 
   
   % print progress 
   progress_text= sprintf('analyzing expt %d of %d \n',index_num, length(idx));
    progress_text_length = length(progress_text);
    fprintf(progress_text);
   
        curr_idx = idx(index_num,:);
        
       
        
        
    for expt = 1: length(Passive.DataDirs)
        
        
        
        % get corr list for  experiment
        try
            Lcorr_expt = Passive.Corrs{1,1}.LCorr(expt,:);
            Ncorr_expt = Passive.Corrs{1,1}.NCorr(expt,:);
        catch 
            Lcorr_expt = Passive.Corr.LCorr(expt,:);
            Ncorr_expt = Passive.Corr.NCorr(expt,:);
        end 
        
        
        %get CellIDs for each index and each experiment
        idx1 = indices{curr_idx(1)}(Passive.experiment_list == expt);
        idx2 = indices{curr_idx(2)}(Passive.experiment_list == expt);
        % convert indices to numbers
        
        
        % subset list of all corrs to just look at corrs between two
        % classes
      
        LCorr(expt,:) = cellfun(@(x)  x(idx1,idx2),Lcorr_expt,...
                        'UniformOutput',0);
        NCorr(expt,:) = cellfun(@(x)  x(idx1,idx2),Ncorr_expt,...
                        'UniformOutput',0);
                    
     
     
      
        
        
    end
    % Statistics and  Significance
    [N_mu,N_std,N_CI,N_sig] = getCorrStatistics(NCorr,'Noise2'); 
    [L_mu,L_std,L_CI,L_sig] = getCorrStatistics(LCorr,'Signal'); 
    
    % Plotting 
     names_text = sprintf('%s, %s',...
                                    names{curr_idx(1)},...
                                    names{curr_idx(2)});
  
   % plot and save signal 
    fig_name_signal =[names_text ': Signal']   ;               
    Bars_by_level(L_mu,L_CI,fig_name_signal,L_sig);
    if ~isempty(savepath) 
        saveCorrFig(fig_name_signal,savepath)
    end
    
  %plot and save noise 
    fig_name_noise =[names_text ': Noise'] ;
    Bars_by_level(N_mu,N_CI, fig_name_noise,N_sig)
    if ~isempty(savepath) 
        saveCorrFig(fig_name_noise,savepath)
    end
    
    % Packaging
    Out = [];
    Out.LCorr = LCorr;
    Out.NCorr = NCorr;
    Out.Signal_Sig = L_sig;
    Out.Noise_Sig  = N_sig;
    Out.Signal_Stats.mean = L_mu;
    Out.Signal_Stats.std = L_std;
    Out.Noise_Stats.mean = N_mu;
    Out.Noise_Stats.std = N_std;
                     
                                
    CorrsBetween{index_num,1} = Out;
    CorrsBetween{index_num,2} = names_text;
    
end

function saveCorrFig(fig_name,savepath)
        fig_name =  matlab.lang.makeValidName(fig_name)
    outpath = fullfile(savepath,[fig_name, '.pdf']);
    print(outpath,'-dpdf','-bestfit')
    pause(.1)
    close(gcf)
