if not exist(Passive,'var')
    load('//vault3/data/kelson/Passive/')
end


names = {'Tone_on','Noise_on','Tone_off'};
DFF_mu = squeeze(nanmean(Passive.DFF,2));
[~,onsets]  = max(DFF_mu);

Noise_idx = onsets > 30 & onsets <= 60;
Tone_idx  = onsets > 60 & onsets <= 90;
Offset_idx = onsets > 90 & onsets<= 120;

indicies = {Tone_idx,Noise_idx,Offset_idx};
 idx = nchoosek( 1:length(indicies),2 );
    
for index_num =  1:length(idx);
   
        curr_idx = idx(index_num,:);
    for expt = 1: length(Passive.DataDirs);
        
        
        
        % get corr list for  experiment
        
        corr_expt = Passive.Corrs{1,1}.Corr{expt};
        
        
        
        
        %get CellIDs for each index and each experiment
        idx1 = indicies{curr_idx(1)}(Passive.experiment_list == expt);
        idx2 = indicies{curr_idx(2)}(Passive.experiment_list == expt);
        % convert indices to numbers
        
        
        
        expt_corrs{expt} = corr_expt(idx1,idx2);
        
    end
    CorrsBetween{index_num,1} = expt_corrs;
    CorrsBetween{index_num,2} = sprintf('%s, %s',...
                                    names{curr_idx(1)},...
                                    names{curr_idx(2)});
    
end

%%%% Option 2 plot all in 3x3
   