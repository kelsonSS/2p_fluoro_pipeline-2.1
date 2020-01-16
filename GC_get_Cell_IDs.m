function CellIDs = GC_get_Cell_IDs(dataDir,Type)


% DataDir - a cell of filepaths to the folders containing the fluorescence,
% experimental, and cell position data for each experiment

% Type determines hwo the data is to be segmented current options
% 'Passive' - No segmentation is performed
% 'Passive-SNR' - segments the data for each level 
% 'Passive-TN'   

  if ~exist('Type','var')
        Type = 'Passive';
    end

dir_main = 'G:\My Drive\2p_fluoro_pipeline 2.1\Granger';
dir_data = '\\vault3\Data\Kelson\Analyzed';

dir_funcs = strcat(dir_main,'Function Directory');
dir_out_main = strcat('G:\My Drive\','GrangerResults\');
% Add directories to path
%addpath(dir_data,dir_funcs)

% Load files from directory if Datadir is empty
if ~exist('dataDir','var')
    gfile_c = 1 ; % gfile counter
    dataDir = {};
    while true
        dataDir{end+1} = uigetdir(dir_data)
        an =  input('Continue? [1/0]');
        if an == 1
            continue
        else
            break
        end
    end
end


%%
irun = 1;
for iif = 1:length(dataDir) %For each file in a directory
    
    % file directory setup
    Datapath = dataDir{iif};
    fn = strrep(Datapath,'-','_');
    fn = strsplit(fn,'\');
   
    
    % grab and Preprocess data
    Data = Fluoro_to_Table(dataDir{iif});
    
    cdef = Data.CellID{1};
    xy = cdef.ptsIdx(:,2:3);
    Fluoro = Data.DFF;
    FLO = Data.FreqLevelOrder;
    
    %%% Data Formatting stage
    %eHDFF = Data.eHDFF; % All Trial Responses for early Hits
    %dHDFF = Data.dHDFF; % All Trial Responses for delayed Hits
    %     Nfreq = numel(Data.PDFF) % Number of freq tones
    %     [N,R,Ncellst] = size(Data.PDFF{1});
    %     PDFF = Data.PDFF{1};
    %     for ff = 2:Nfreq
    %         PDFF = cat(2,PDFF,Data.PDFF{ff});
    %     end
    %     size(PDFF)
    %     PDFF = zeros(N,Nfreq*R,Ncellst);
    %     for ff = 1:Nfreq
    %         PDFF(:,R*(ff-1)+1:R*ff,:) = Data.PDFF{ff};
    %     end
    
    %
    %     figure
    %     plot(reshape(mean(PDFF,2),[size(PDFF,1),size(PDFF,3)]))
    %
    % Choose different modes e.g. 'Hit' for Hits and 'Miss' for Miss
    % Currently must comment out which version you are not using
    
    % TODO- if code is functionalized, make this a case-switch
  
    switch Type
        case 'Passive'
            % Passive-(ALL Data)
            mode = {'Passive'}; nmode = numel(mode);
            resp{1} = Fluoro;
            ID = [1:size(Fluoro,3)];
            
        case 'SNR'
            % Passive SNR
            mode = {'Tone','+30db SNR','+20db SNR','+10db SNR'};
            nmode = numel(mode);
            Levels =table2array(Data.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
            for lvl = 1:length(uL)
                  resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
            end 
            ID = [1:size(Fluoro,3)];
            
        case 'Tones'
            % Passive-SNR- Tone Responsive Neurons
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
            Tone_idx  = onsets > 60 & onsets <= 90;
            Fluoro = Fluoro(:,:,Tone_idx);
            ID = find (Tone_idx);
            
            mode = {'TR- Tone','TR +30db SNR','TR +20db SNR','TR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(Data.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                  resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end 
            
        case 'Noise'
            % Passive SNR Noise Responsive Neurons
            
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
            Noise_idx  = onsets > 30 & onsets <= 60;
            Fluoro = Fluoro(:,:,Noise_idx);
            ID = find(Noise_idx);
            
            mode = {'NR- Tone','NR +30db SNR','NR +20db SNR','NR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(Data.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                 resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end
              
        case 'Offset'   
             % Passive SNR Noise Responsive Neurons
            
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
            Offset_idx  = onsets > 90 & onsets <= 120;
            Fluoro = Fluoro(:,:,Offset_idx);
            ID = find(Offset_idx);
            
            mode = {'NR- Tone','NR +30db SNR','NR +20db SNR','NR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(Data.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                 resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end
              
              
    end
    
    
    
    
   
   
   
   %%
 %%  
 
    
    
    %clear Data 
    clear Fluoro
    clear Data
    
    
    %%% Locate and remove Nan samples from recordings
    ccorr = [];
    for imd = 1:nmode
        imd
        Resp = resp{imd};
        Rnan = isnan(Resp);
        cnan = mean(mean(Rnan,1),2);
        ccorr(imd,:) = squeeze(cnan == 1);
        %inan = sum(sum(Rnan,1),3);
        %rcorr = find(inan)
        %resp{imd}(:,rcorr,:) = [];
        %RR(imd) = RR(imd) - numel(rcorr)
    end
    cpass = find(sum(ccorr,1) == 0);
    Ncellst = numel(cpass);
    Ncellz = 20;
    Ncells = min(Ncellst,Ncellz);
    
    for imd = 1:nmode
        Resp = resp{imd}(:,:,cpass);
         Rnan = isnan(Resp);
         inan = sum(sum(Rnan,1),3);
         rcorr = find(inan == 0);
         respc{imd} = Resp(:,rcorr,:);
    end
    xy = xy(cpass,:);
    
    %%% Cell Selection need to 10 - 20 cells, depending on number of
    %%% trials
    % Select a Subset of Cells with highest response variability
    RR = zeros(nmode,1);
    varR = zeros(Ncellst,nmode); 
    
    Grand = 1;  
    
    if  Ncellst > Ncellz 
           
        for imd = 1:nmode
            Resp = respc{imd};
            RR(imd) = size(Resp,2); 
            for c = 1:Ncellst
                varR(c,imd) = mean(var(Resp(:,:,c)));
            end        
        end
       % The total variance, varS , for each cell is equal to the weighted average
       %  of variances of a given cell across conditions
        varS = (varR*RR)/sum(RR);
        %varS2 = (varR(:,1:2)*RR(1:2))/sum(RR);

        [~,indm] = sort(varS,'descend');    
        cellids = sort(indm(1:Ncells))  ;   
    else
       cellids = [1:Ncellst];
       
    end
    
    CellIDs{iif} = ID(cellids);
end 