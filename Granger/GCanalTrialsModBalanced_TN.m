%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% GC Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Results = GCanalTrialsModBalanced_TN(TN,Type,RedChannel)

% TN - an experiment object created by Fluoro_to_table
% that conrains the cell of filepaths to the folders containing the fluorescence,
% experimental, and cell position data for each experiment

% Type determines hwo the data is to be segmented current options
% 'Passive' - No segmentation is performed
% 'Passive-SNR' - segments the data for each level 
% 'Passive-TN'   

% number of neurons to include = (.1 * Timepoints * repeats )/ time-windows

% Red channel is a boolean that, if true, restricts the data to only cells 
% that were labeled in the Data.RedCellID
tic;
  if ~exist('Type','var')
        Type = 'SNR'
  end
  
  if isnumeric(Type)
      % the type is stored as the first number in the array and the rest of
      % the numbers are the index of which cells are of that given type. We
      % store and seperate that out here. 
    ClassNum = Type(1);
    TypeName = TN.Classes{Type(1)};
    Type = Type(2:end);
  else 
      TypeName = Type;
  end 
  
  
  if ~exist('RedChannel','var')
      RedChannel = false; 
  end 
  

%  dir_main = 'G:\My Drive\2p_fluoro_pipeline 2.1\Granger';
% dir_data = '\\vault3\Data\Kelson\Analyzed';
% 
% dir_funcs = strcat(dir_main,'Function Directory');
dir_out_main = strcat('D:\Analyzed\GrangerResults\');
% Add directories to path
%addpath(dir_data,dir_funcs)


try
    dataDir = TN.DataDirs;
catch
    try
        dataDir =TN;
        clear TN
    catch
    end 
end 

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
for expt = 1:length(dataDir) %For each file in a directory
    
%     % file directory setup
%     try
%     Datapath = dataDir{iif};
%     catch 
%         Datapath = dataDir;
%     end
   
 %% grab and Preprocess data
    try
     %% data unpacking
     Data = dataDir{expt};
     cdef = TN.CellID{expt};
     xy = cdef.ptsIdx(:,2:3);
     
     if ~isstruct(Data) && ischar(Data)
         Datapath = Data;
         fn = strrep(Datapath,'-','_');
         fn = strrep(Datapath,'/','\');
         fn = strsplit(fn,'\');
         fileName1 = strcat(fn{end-2} ,'_',fn{end-1}, '_' , TypeName);
         if RedChannel
             fileName1 = [fileName1 '_Labeled'];
         end  
         
         
         dir_out = fullfile(dir_out_main,fileName1);
         
             
         dir_out_all{expt} = dir_out;
         
     
         
     
     if ~isempty(dir(fullfile(dir_out,'DataResults*')))
         continue
     end     

      
 
  
     end 
     
  
     
     %% name and file path creation
  
     fn = strrep(Data,'-','_');
     fn = strrep(Data,'/','\');
     fn = strsplit(fn,'\');
     fileName1 = strcat(fn{end-2} ,'_',fn{end-1}, '_' , TypeName);
     
         
     dir_out_all{expt} = dir_out;    
     
     if ~isempty(dir(fullfile(dir_out,'DataResults*')))
         continue
     end 
     
     
     if ~exist(dir_out,'dir')
         mkdir(dir_out)
     end 
 
   
     
    catch
        % if given empty cell
        continue
    end 
     
    try
        expt_idx = TN.experiment_list == expt;
        clean_idx = TN.Clean_idx(expt_idx);
        active_idx = TN.active{expt_idx,2}>0;
        Fluoro = TN.DFF_Z(:,:, expt_idx);
        Fluoro = Fluoro(:,:,clean_idx & active_idx);
        
       
    catch
        warning('could not use normal variables, mmanual inspeciton needed')
    
        Data = Fluoro_to_Table(Data);   
        FLO = Data.FreqLevelOrder;
           
        cdef = Data.CellID{1};
        xy = cdef.ptsIdx(:,2:3);
        active_idx = Data.active{:,2} >= 1; % find all cells active P<.05
     try
        Fluoro = Data.DFF(:,:,active_idx);
      
    
    catch 
        Fluoro = Data.DFF{1}(:,:,active_idx);  
    end 
 
    
    DFF_Z = squeeze( ( Fluoro - nanmean(nanmean(Fluoro )) ./ nanstd(nanstd(Fluoro))));
    
    % remove artifactual data by testing if the silence frames have
    % unusually high variance 
 
    
    clean_idx = squeeze(max(abs(nanstd(DFF_Z(1:30,:,:),[],2)))<20 ); 
    
    Fluoro = Fluoro(:,:,clean_idx);
    
    end
    
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

  
    if  isnumeric(Type)
               % Use predetermined Clusters  
               
            DFF_mu = squeeze(nanmean(Fluoro,2));
            idx =  TN.Class_idx == ClassNum ; 
            idx = idx(TN.experiment_list == expt);
            idx = idx(clean_idx & active_idx);
            
            try
                Fluoro = Fluoro(:,:,idx);
            catch 
                Fluoro = Fluoro(:,:,idx(1:size(Fluoro,3)) );
            end 
            
            
            
            mode = { ' Tone',' +30db SNR',' +20db SNR',' +10db SNR'};
            mode = cellfun(@(x) [TypeName ' '  x],mode,'UniformOutput' ,0);  
    
            
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                 resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end    
        
        
        
    else 
    switch Type
        
      
        case 'Passive' 
            % Passive-(ALL Data)
            mode = {'Passive'};
            nmode = numel(mode);
            resp{1} = Fluoro;
            
        case 'Active'
            % Passive-(ALL Data)
            mode = {'Passive'};
            nmode = numel(mode);
            resp{1} = Fluoro;
  
        case {'Cluster_1','Cluster_2','Cluster_3','Cluster_4'}
               % Passive SNR
            clust_num = str2num(Type(end));
            num_neurons = sum(clusters == clust_num); 
            fprintf('Experiment %d,Cluster %d: %d neurons \n'...
                    ,expt,clust_num,num_neurons)
                if num_neurons <10
                    continue
                end 
            
            freqs = unique(FLO(:,2));
            mode = cellfun(@num2str , table2cell(freqs),'UniformOutput',0);
           % mode = {'Tone','+30db SNR','+20db SNR','+10db SNR'};
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
            
            clear resp
            for lvl = 1:length(uL) 
                  mode{lvl} = sprintf('Cluster_%s_%d_db', TN.Classes{clust_num},uL(lvl));
                  resp{lvl} = Fluoro(:,Levels == uL(lvl),clusters == clust_num);
                 
            end 
             
            
            
        case 'SNR'
            % Passive SNR
            freqs = unique(FLO(:,2));
            mode = cellfun(@num2str , table2cell(freqs),'UniformOutput',0);
           % mode = {'Tone','+30db SNR','+20db SNR','+10db SNR'};
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
            for lvl = 1:length(uL)
                  resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
            end 
            
           
        case 'EachFreqLevel'
            % collect and seperate unique freq/level indicies 
            Levels =table2array( TN.FreqLevelOrder(:,2) );
             Freqs =table2array( TN.FreqLevelOrder(:,1) );
              uL = unique( Levels );
              uF = unique( Freqs );
              % iterate over them to create groups 
              mode ={};
              resp_idx = 1;
             for freq = 1:length(uF)
                 for lvl = 1:length(uL)
                     mode{end+1} = sprintf('%d hz, %d dB Tone',uF(freq),...
                                            uL(lvl) );
                     Fluoro_idx = Levels == uL(lvl) & Freqs == uF(freq);
                                        
                     resp{resp_idx} = Fluoro(:,Fluoro_idx,:);      
                     
                     resp_idx = resp_idx + 1; 
                                        
                     
            
                 end
             end 
             nmode = resp_idx - 1; 
   
            
        case 'HitMissEarly'
            mode = {'Hit','Miss','Early'}
            nmode = numel(mode);
            
            hit_idx = logical(TN.handles.Hits);
            early_idx = logical(TN.handles.Early);
            miss_idx = logical(TN.handles.Miss);
            
            resp{1} = Fluoro(:,hit_idx,:);
            resp{2} = Fluoro(:,early_idx,:);
            resp{3} = Fluoro(:,miss_idx,:);
         
            
       case 'HitMiss'
            mode = {'Hit','Miss'}
            nmode = numel(mode);
            
            hit_idx = logical(TN.handles.Hits);
            early_idx = logical(TN.handles.Early);
            miss_idx = logical(TN.handles.Miss) | logical(early_idx) ;
            assert(~any(miss_idx>1))
            
            resp{1} = Fluoro(:,hit_idx,:);
        
            resp{2} = Fluoro(:,miss_idx,:);
            
        case 'HitLevels'
            levels = TN.handles.Levels; 
            uLevels = TN.handles.uLevels;
             
             mode = arrayfun(@num2str,uLevels,'UniformOutput',0);
             nmode = numel(mode);
             hit_idx = logical(TN.handles.Hits);
            
             for lvl = 1:length(uLevels)
                 
                 lvl_idx = levels ==uLevels(lvl);
                 main_idx = lvl_idx' & hit_idx;
                 
                 resp{lvl} = Fluoro(:,main_idx,:);
             end 
            
        case 'Tones'
            % Passive-SNR- Tone Responsive Neurons
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
            Tone_idx  = onsets > 60 & onsets <= 90;
            Fluoro = Fluoro(:,:,Tone_idx);
            
            mode = {'TR- Tone','TR +30db SNR','TR +20db SNR','TR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
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
            
            mode = {'NR- Tone','NR +30db SNR','NR +20db SNR','NR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                 resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end
              
        case 'Offset'   
             % Passive SNR Tone-OFF/Ramping Responsive Neurons
            
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
            Offset_idx  = onsets > 90 & onsets <= 120;
            Fluoro = Fluoro(:,:,Offset_idx);
            
            mode = {'NR- Tone','NR +30db SNR','NR +20db SNR','NR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                 resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end
       case 'Off'   
             % Passive SNR Noise-OFF Responsive Neurons
            
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
             Off_idx  = onsets > 120;
            Fluoro = Fluoro(:,:,Off_idx);
            
            mode = {'NR- Tone','NR +30db SNR','NR +20db SNR','NR +10db SNR'};
            nmode = numel(mode);
            Levels =table2array(TN.FreqLevelOrder(:,2));
            uL = unique(Levels);
            
              for lvl = 1:length(uL)
                 resp{lvl} = Fluoro(:,Levels == uL(lvl)  ,:);
              end    
              
              
    end
    end 
    
    
    
   
   
   
   %%
 %%  
 
    
    
%     %clear Data 
%     clear Fluoro
%     clear Data
%     
%     

  if ~exist('clust_num','var')
      clust_num = 1;
  end 
    
  
    ccorr = [];
       %%% Locate and remove Nan samples from recordings
    cpass = ones(size(resp{1},3),1);
       for imd =1:nmode
        
        mode_pass =  ~squeeze(any(any(isnan(resp{imd}))));
           if sum(mode_pass) == 0;
               resp_temp = resp{imd};
               resp_temp(isnan(resp_temp)) = 0;
               resp{imd} = resp_temp;
               clear resp_temp
                 mode_pass =  ~squeeze(any(any(isnan(resp{imd}))));
           end 
%         Resp = resp{imd};
%         Rnan = isnan(Resp);
%         cnan = mean(mean(Rnan,1),2);
%         ccorr = cat(1,ccorr, squeeze(cnan == 1));
        %inan = sum(sum(Rnan,1),3);
        %rcorr = find(inan)
        %resp{imd}(:,rcorr,:) = [];
        %RR(imd) = RR(imd) - numel(rcorr)
        cpass = cpass & mode_pass;
        
       
    end
    %cpass = find(sum(ccorr,1) == 0);
   Ncellst = sum(cpass);
    % Ncellst = numel(cpass);
   
        disp(num2str(Ncellst))
  
     %Ncells = min(Ncellst,Ncellz);
  % Test to see if responses are okay with using all active cells
     %
    % Ncellz = sum(active_idx) 
     Ncellz = 20;
     Ncells = min(Ncellst,Ncellz);
     
    for imd =1:nmode
        Resp = resp{imd}(:,:,cpass);
         Rnan = isnan(Resp);
         inan = sum(sum(Rnan,1),3);
         rcorr = find(inan == 0);
         respc{imd} = Resp(:,rcorr,:);
    end
    try
    xy = xy(cpass,:);
    catch
        continue
    end 
    %%% Cell Selection need to 10 - 20 cells, depending on number of
    %%% trials
    % Select a Subset of Cells with highest response variability
    RR = zeros(nmode,1);
    varR = zeros(Ncellst,nmode); 
    
    Grand = 1;  
    
    if  Ncellst > Ncellz 
           
        for imd =1:nmode
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
        cellids = [1:Ncells]';
        RR =  cellfun( @(x) size(x,2), respc);
        if length(cellids)<10
            continue
        end
    end    
    
   save(fullfile(dir_out,['CellIDs']),'cellids') 
    %%% Random Trial Selection [for balanced analysis]
    % If two different modes, balance the number of trials. e.g. hit vs
    % miss
    indc = find(RR < 10) ;
     if numel(indc) == 0
         [Rmin,imin] = min(RR);
         respx = respc;
         icc = 1:nmode; icc(imin) = [];
         for i = icc
             respx{i} = respc{i}(:,sort(randperm(RR(i),Rmin)),:);
         end

      % Uncomment here if u want to run random trial selection  
%     elseif numel(indc) == 1
%         [RRm,im] = sort(RR,'ascend')
%         icc = 1:nmode; icc(im(1)) = [];
%         mode = mode(icc) 
%         nmode = numel(mode); 
%         respx = respc(icc); 
%         for i = 1:nmode
%             respx{i} = respx{i}(:,sort(randperm(size(respx{i},2),RRm(2))),:);
%         end
     else
         % if min # of trials is less than 10, GC analysis is not reliable
%         continue
     end


%% init variables 
    Devtotal = zeros(Ncells,Ncells,nmode); DLLtotal = Devtotal;
    wftotal = cell(Ncells,nmode);
    robstotal = cell(Ncells,nmode);
    sig2ftotal = zeros(Ncells,nmode);
    sig2rtotal = zeros(Ncells,Ncells,nmode);
    Bftotal = zeros(Ncells,nmode);
    Brtotal = zeros(Ncells,Ncells,nmode);
    rhatftotal = cell(Ncells,nmode);
    rhatrtotal = cell(Ncells,Ncells,nmode);
    gammaOpttotal = zeros(Ncells,nmode);
    cellidstotal = cell(nmode,1);
 
    %% main loop 

    
    
    for imd =1:nmode
    fprintf('Expt %d/%d mode:%d/%d \n', expt,length(dataDir),imd,nmode)
    Resp = respx{imd}(:,:,cellids);
    [N,R,Ncells] = size(Resp);
    %%% GC Analysis: Multi-trial Multi-Cell Mod
    % Estimation Setting:
    LL = 1000;
    LR = 10;
    % Cross-History covariate setting
    WH = 5; %window: average the values over the window over moving window: filters data over 6 Hz that are not calcium signals
    Mhc = 5; %nunber of samples per cell in the model
    
    Whc = [1 WH*ones(1,Mhc-1)]; % Cross-Hist Kernel-make sure that the length latency of interaction + latency of Ca maging
    Lhc = sum(Whc);
    
    % Self-Hist covariate setting
    Mhcs = 5;
    Whcs = [1 WH*ones(1,Mhcs-1)]; % Self-Hist Kernel
    Lhcs = sum(Whcs);
    Mf = (Ncells-1)*Mhc+Mhcs; Mr = (Ncells-2)*Mhc+Mhcs;
    %%% Cross-Validation Step
    gammacv = [0:0.2:3]; % check for 15 diff gamma for each cell
    Whcv = {Whc}; Whscv = {Whcs};
    LLcv = 100; %LL
    LRcv = LR;
    [gammaOpt,Jcost] = CrossValidModTrial2(Resp,gammacv,Whcv,Whscv,LLcv,LRcv);
    %save(['CrossValidData' fileName 'mode' mode{imd} 'Ncells35Wh' num2str(WH) 'Md' num2str(Mhc) 'ModWhc2.mat'],'gammaOpt','gammacv','Whcv','LLcv','LRcv')
    gammaOpttotal(:,imd) = gammaOpt;
    Lm = max([Lhc,Lhcs]);
    Np = N - Lm;
    %%% Form Cross-History Covariate Matrices
    Xcz = zeros(Np,Mhc,Ncells,R);
    for r = 1:R
        Xcz(:,:,:,r) = FormHistMatrixv2(squeeze(Resp(:,r,:)),Whc,Lm); % resp of size [N x Ncells]   
    end
    %%% Form Observation arrays for all cells
    Robs = Resp(Lm+1:end,:,:); %Robs = Resp; 
    % Zero-mean the observation vectors
    % This will remove the need for intercept (mu) estimation
    for cc = 1:Ncells
        Robs(:,:,cc) = Robs(:,:,cc) - ones(Np,1)*mean(Robs(:,:,cc));
    end
    for ct = 1:Ncells % ct: index of target cell
        N_indto = num2str(ct);
        cellC = 1:Ncells; cellC(ct) = [];
        gamma = gammaOpt(ct);
        %%% Form SELF-History Covriate matrix
        Xcsz = zeros(Np,Mhc,R);
       for r = 1:R
            Xcsz(:,:,r) = FormHistMatrixv2(squeeze(Resp(:,r,ct)),Whcs,Lm); % resp of size [N x Ncells]   
        end
        %%% Form Effective response vector
        robsz = Robs(:,:,ct);
        reff = robsz(:); Neff = numel(reff);
        %reffn = reff/sqrt(var(reff));
        %%% Full Model
        %%% Form standardized Full Model Covariate Matrix
        [Xneff,Dnf] = FormFullDesignv3(Xcz,Xcsz,ct);
        %%% Filtering/Estimation Setting
        tic
        [wnf,sig2f,Devf,Bf] = StaticEstim(Xneff,reff,gamma,LL,LR);
     
        %%% Save estimated parameters/arrays
        % Save Tuning AR Coeffs full model for each cell
        wftotal{ct,imd} = wnf;
        sig2ftotal(ct,imd) = sig2f;    
        robstotal{ct,imd} = robsz;
        %Dtotal{ct} = Dnf; % SAVE scalings for cells
        % Reconstructed observation vector full model
        rhatf = zeros(Np,R);
        for r = 1:R
            Xnf = Xneff((r-1)*Np+1:r*Np,:);
            rhatf(:,r) = Xnf*wnf;
        end
        rhatftotal{ct,imd} = rhatf;
        Bftotal(ct,imd) = Bf;
        
        for cf = cellC % cf : index of Electrode/Unit we check causality from to the target 'ct'
            N_indfrom = num2str(cf) ;
            
            
            line_length =  fprintf(['Estimating G-Causality from cell ' N_indfrom ' to cell ' N_indto ' ...']);
            %%% Reduced Model
            Xnefr = Xneff;
            cc = find(cellC == cf);
            Xnefr(:,(cc-1)*Mhc+Mhcs+1:cc*Mhc+Mhcs) = [];
            %%% Estimation Setting           
            [wnr,sig2r,Devr,Br] = StaticEstim(Xnefr,reff,gamma,LL,LR);
        
            % Save estimated params/arrays
            sig2rtotal(ct,cf,imd) = sig2r;
            % Reconstructed observation vector reduced model
            rhatr = zeros(Np,R);
            for r = 1:R
                Xnr = Xnefr((r-1)*Np+1:r*Np,:);
                rhatr(:,r) = Xnr*wnr;
            end
            rhatrtotal{ct,cf,imd} = rhatr;
            Brtotal(ct,cf,imd) = Br;
            %%% Granger Causality metric: Compute Difference of Deviance Statistic
            Devd = Devf - Devr;
            DB = Bf - Br;
            DLL = Devd - DB;
            Devtotal(ct,cf,imd) = Devd;
            DLLtotal(ct,cf,imd) = DLL;
            fprintf(repmat('\b',1,line_length));
        end
    end
    cellidstotal{imd} = cellids;
end

%%% Statistical Test: J-statistics based on FDR
Md = Mhc; %bigger value test higher stat val.
alpha = 0.01; %sig. level; 

DkSmth = Devtotal -Md;
%[GCJL,GCJH] = FDRcontrolBHv4(Devtotal,DkSmth,alpha,Md,wftotal,Whc);
Mlow = 2 %determines which segemetns of Whx are used to dtermine inhibitory vs ecit: limit length of window to avoid phase ambiguoity
[GCJL,GCJH] = FDRcontrolBHv5(Devtotal,DkSmth,alpha,Md,wftotal,Whc,Mlow);

%% Plot Results and Figures

hf = figure;
set(hf, 'position',[10 10 810 610]); 

% find optimal configuration of subplots for most square config 
%( ex ifnmode = 8, rc = [3 3] 
rc = numSubplots(nmode);

for imd =1:nmode
    cgtxt = {cellidstotal{imd}};
    subplot(rc(1),rc(2),imd), imagesc(Devtotal(:,:,imd))
    cmax = max(max(Devtotal(:,:,imd)));
    cmin = 0;
    
    caxis([cmin cmax])
    title( mode{imd} )
   % title([ mode{imd} ' @ Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd))) ])
    set(gca,'XTick',[1:Ncells],'YTick',[1:Ncells]);
    set(gca,'xticklabel',cgtxt,'yticklabel', cgtxt);%{cellid}
    %xtickangle(90)
    colorbar;
end
%title([ mode{imd} ' @ Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd))) ])
%suptitle(['Estimated deviance @ Md = ' num2str(Md) ])
colormap jet
%saveas(gcf,fullfile(dir_out,['FigDevmap' num2str(Md) 'Whc' num2str(WH) 'gammaAvg' num2str(mean(gammaOpt)) 'Modnnrm.fig']));
saveas(gcf,fullfile(dir_out,['FigDevmap' num2str(Md) 'Whc' num2str(WH)  'Modnnrm.fig']));
close(hf)

cmax = 1; cmin = -1;
hf = figure;
set(hf, 'position',[10 10 810 610]); 
for imd =1:nmode
    cgtxt = {cellidstotal{imd}};
    subplot(rc(1),rc(2),imd), imagesc(GCJL(:,:,imd))
    caxis([cmin cmax])
    title( mode{imd} )
    %title([ mode{imd} ' @ FDR =' num2str(alpha) ' Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd)))])
    set(gca,'XTick',[1:Ncells],'YTick',[1:Ncells]);
    set(gca,'xticklabel',cgtxt,'yticklabel', cgtxt);%{cellid}
    %xtickangle(90)
    colorbar;
end
%title([ mode{imd} ' @ FDR =' num2str(alpha) ' Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd)))])
%suptitle(['Detected GC maps @ FDR = ' num2str(alpha) ' Md = ' num2str(Md) ])
colormap jet
saveas(gcf,fullfile(dir_out,['FigGCmapsLowLat_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'ModnnrmRun' num2str(irun) '.fig']));
close(hf)



%%% Plot the GC network maps based on Original location of cells
ymax = ceil(max(max(xy))/100)*100 + 10;
%Set min/max for marker size and thickness 
thick_max = 7;
Marker_max = 10;%Marker_min= 10;
aa = 5;
crds = xy;

hf = figure;
set(hf, 'position',[0 0 800 800]); % To set position and scaling of a figure window
for imd =1:nmode
    XY = xy(cellidstotal{imd},:);
    crdsf = XY; % Set Locations of Selected Cells
    rlen = zeros(Ncells); rtheta = rlen;
    for i = 1:Ncells
        ji = 1:Ncells; ji(i) = [];
        for j = ji
            rlen(i,j) = norm(crdsf(i,:) - crdsf(j,:));
            rtheta(i,j) = atan((crdsf(i,2) - crdsf(j,2))/(crdsf(i,1) - crdsf(j,1)));%abs(crdsf(i,2) - crdsf(j,2))/rlen(i,j);
        end
    end

    hs = subplot(rc(1),rc(2),imd);
    axis square   
    % Plot Markers for All cells first    
    for i = 1:size(crds,1) %v = Channels
        plot(crds(i,1),crds(i,2),'o','MarkerSize',Marker_max,'MarkerFaceColor','w','MarkerEdgeColor','r')
        text(crds(i,1)+5,crds(i,2)+5,num2str(i),'color','r','Fontsize',12)
        hold on
    end
    % Plot Markers for Selec cells    
    pp = Marker_max;
    for i = 1:Ncells %v = Channels
        plot(crdsf(i,1),crdsf(i,2),'o','MarkerSize',pp,'MarkerFaceColor','y','MarkerEdgeColor','r')
        %text(crdsf(i,1)+10,crdsf(i,2)+10,num2str(Channels(i)),'Fontsize',10)
    end
    % Plot GC links
    for r = 1:Ncells
        for c = 1:Ncells
            if c ~= r
                ll= abs(GCJL(r,c,imd))*thick_max ;%+ thick_min
                if ll> thick_max, ll = thick_max; end
                if (ll>0)
                    if r > c, aar = -aa;
                    else aar = aa;
                    end
                    if sign(GCJL(r,c,imd))< 0
                        cll = 'b';
                    else
                        cll = 'r';
                    end
                    offs = aar*[-sin(rtheta(r,c)) cos(rtheta(r,c))];
                    arrow([crdsf(c,1)+offs(1) crdsf(c,2)+offs(2)],[crdsf(r,1)+offs(1) crdsf(r,2)+offs(2)],'Color',cll,'Width',ll,'Length',15,'TipAngle',30);
                end
            end
        end
    end
    axis([0 ymax 0 ymax])
    axis xy
    set(gca,'xticklabel',{[]}, 'yticklabel', {[]});
    ax = gca;
    ax.Color = 'k';
    title(['GC net maps @ mode ' mode{imd}])
end
saveas(gcf,fullfile(dir_out,['FigGCnetworkMaps_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'ModCVModWh2ModnnrmRun' num2str(irun) '.fig']));
close(hf)

    %%% Compute different Statistic measures
    GCnumbers = zeros(nmode,1);
    GClengths = cell(nmode,1);
    GCangles = cell(nmode,1);
    GCavgL = zeros(nmode,1); %GCavgLw = GCavgL;
    
    for imd =1:nmode
        %%% Compute number of GC links
        GCnumbers(imd) = numel(find(GCJL(:,:,imd)));
        %%% Compute GC link lengths
        GCL = (GCJL(:,:,imd) ~= 0).*rlen;
        GClengths{imd} = GCL(find(GCL));
        GCavgL(imd) = mean(GCL(find(GCL)));
        %GCLw = GCJL(:,:,imd).*rlen;
        %GCavgLw(imd) = mean(GCLw(find(GCLw)));
        %%% Compute GC link angles
        [ii,jj] = find(GCJL(:,:,imd));
        GCtheta = zeros(numel(ii),1);
        for l = 1:numel(ii)
            dy = xy(cellids(ii(l)),2) - xy(cellids(jj(l)),2);
            if dy*rtheta(ii(l),jj(l)) < 0
                GCtheta(l) = rtheta(ii(l),jj(l)) - sign(rtheta(ii(l),jj(l)))*pi;
            else
                GCtheta(l) = rtheta(ii(l),jj(l));
            end
        end
        GCangles{imd} = GCtheta;
    end
    save(fullfile(dir_out,['DataStatsResultsGCanal_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'Modnnrm.mat']),'GCnumbers','GClengths','GCangles','GCavgL','mode') %,'GCavgLw'   

    %%% Plot activation/inhibition effects
    % Plot responses associated with two target cells: Confirm detected GCs visually
    nsubc = 6;
    for imd =1:nmode
        nGCs = numel(find(GCJL(:,:,imd)))
        nsubr = ceil(nGCs/nsubc); %min(cnt,nGCs)
        
        [Ncells,~,~] = size(Devtotal);
        %cids = cellidstotal{imd};
        Dev = Devtotal(:,:,imd); Devv = Dev(:);
        [~,indev] = sort(Devv,'descend');
        ii = zeros(nGCs,1); jj = ii;
        for l = 1:nGCs
            jj(l) = ceil(indev(l)/Ncells);
            ii(l) = indev(l) - (jj(l)-1)*Ncells;
        end
        GClist = [ii,jj]
        
        ymax = max(max(mean(respx{imd}(:,:,cellids),2)));
        ymin = min(min(mean(respx{imd}(:,:,cellids),2)));        
        hfimd = figure;
        for l = 1:nGCs
            ct = GClist(l,1);
            cf = GClist(l,2);
            cellids = cellidstotal{imd};
            respcf = respx{imd}(:,:,cellids(cf)); respct = respx{imd}(:,:,cellids(ct));
            [N,~] = size(respct);
            subplot(nsubr,min(nsubc,nGCs),l), plot(mean(respct,2))
            hold on
            plot(mean(respcf,2),'r')
            title(['GC(' num2str(cellids(cf)) '\rightarrow' num2str(cellids(ct)) ')'])
            xlim([1,N])
            ylim([ymin ymax])
        end
        saveas(hfimd,fullfile(dir_out,['FigEffectResp_' fileName1 mode{imd} 'mode_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'ModErlyDlyModnnrmRun' num2str(irun) '.fig']));
        close(hfimd)
    end

    save(fullfile(dir_out,['DataResultsGCanal_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'ModCVModWh2Modnnrm.mat']),'Devtotal','GCJL','GCJH','wftotal','Whc','Whcs','sig2ftotal','sig2rtotal','robstotal','rhatftotal','rhatrtotal','Bftotal','Brtotal','cellidstotal','gammaOpttotal','LL','LR','alpha')

end 


%% Gather GC data Back into Experiment object 
    Results = [];
    try
    Results = Collect_GC_Results(dir_out_all);
    catch
       warning('Error in collecting GC results!')
     end 
%% [Not neccessary to run!] More Detailed Tests on different parts of analysis

%% Check estimated hist kernels of Self + Cross couplings
% imd = 1
% ct = 11
% 
% wf = wftotal{ct,imd};
% 
% figure
% stem(wf)
% 
% numel(find(abs(wf) > std(wf)/5))


%%
% ct = 11
% 
% rct = robstotal{ct,imd};
% 
% figure
% plot(rct)

%% Find Most Significant GC Links
% imd = 1;
% [Ncells,~,~] = size(Devtotal);
% cids = cellidstotal{imd};
% Dev = Devtotal(:,:,imd);
% Devv = Dev(:);
% [~,indev] = sort(Devv,'descend');
% 
% cnt = 20;
% ii = zeros(cnt,1); jj = ii;
% for l = 1:cnt;
%     jj(l) = ceil(indev(l)/Ncells);
%     ii(l) = indev(l) - (jj(l)-1)*Ncells;
% end
% 
% GClist = [ii,jj]
% GClistp = [cids(ii),cids(jj)]

%% Plot responses associated with two target cells: Confirm detected GCs visually
% imd = 1
% ct = 14
% cf = 5
% 
% respcf = robstotal{cf,imd};
% respct = robstotal{ct,imd};
% [Np,~] = size(rhatftotal{1,1});
% 
% figure
% plot(mean(respct,2)/sqrt(var(respct(:))))
% hold on
% plot(mean(respcf,2)/sqrt(var(respct(:))),'r')
% plot(1:Np,mean(rhatftotal{ct,imd},2),'k','linewidth',2)
% plot(1:Np,mean(rhatrtotal{ct,cf,imd},2),'g','linewidth',2)
% title(['Mean activity of cells ' num2str(cf) ' & ' num2str(ct) ' across trials'])
% legend(num2str(ct),num2str(cf))
% xlim([1,Np])
%saveas(gcf,['Sim1SignifGCs_MeanCellAct' num2str(cf) '&' num2str(ct) '.fig']);
%close(gcf)

% iif = 37
% fileName = files(fileIndex(iif)).name
% load(fileName)
%HDFF = Data.eHDFF; % Hit Trial Responses
%MDFF = Data.MDFF; % Miss Trial Responses
%xy = Data.xy; % contain locations of cells
%respx{1} = HDFF; respx{2} = MDFF;
%clear Data HDFF MDFF

% cellids = cellidstotal{imd};
% respcf = respx{imd}(:,:,cellids(cf)); respct = respx{imd}(:,:,cellids(ct));
% [N,~] = size(respct);
% figure
% plot(mean(respct,2))
% hold on
% plot(mean(respcf,2),'r')
% title(['Mean activity of cells ' num2str(cf) ' & ' num2str(ct) ' across trials'])
% legend(num2str(ct),num2str(cf))
% xlim([1,N])
% 
% figure
% plot(respct,'b')
% hold on
% plot(respcf,'r')
% title(['Activity of cells ' num2str(cf) ' & ' num2str(ct) ' within all trials'])
% xlim([1,N])
% saveas(gcf,['Sim1SignifGCs_CellAct' num2str(cf) '&' num2str(ct) '.fig']);
%close(gcf)

%%
% imd = 2
% ct = 10
% cf = 11
% 
% rhatct = rhatftotal{ct,imd};
% rhatcf = rhatftotal{cf,imd};
% 
% figure
% subplot(1,2,1),plot(rhatct,'b'), hold on, plot(mean(rhatct,2),'k','linewidth',2)
% title(['Activity across trials @' mode{imd} ' ct=' num2str(cellids(ct))])
% xlim([0 N-Lm])
% subplot(1,2,2),plot(rhatcf,'r'), hold on, plot(mean(rhatcf,2),'k','linewidth',2)
% title(['Activity across trials @' mode{imd} ' cf=' num2str(cellids(cf))])
% xlim([0 N-Lm])
% 
% rr = 9
% figure
% subplot(1,2,1),plot(rhatct(:,rr),'b'), hold on, plot(mean(rhatct,2),'k','linewidth',2)
% title(['Activity @' mode{imd} ' ct=' num2str(cellids(ct))])
% xlim([0 N-Lm])
% subplot(1,2,2),plot(rhatcf(:,rr),'r'), hold on, plot(mean(rhatcf,2),'k','linewidth',2)
% title(['Activity @' mode{imd} ' cf=' num2str(cellids(cf))])
% xlim([0 N-Lm])
% 
% %%
% ct = 16
% wnf = wftotal{ct,imd};
% 
% figure
% stem(wnf)
% 
% %% Plot the sample figure of activations 
% hf = figure;
% l = [1,2,4,8]
% for i = 1:4
%     ct = GClist(l(i),1);
%     cf = GClist(l(i),2);
%     cellids = cellidstotal{imd};
%     respcf = respx{imd}(:,:,cellids(cf)); respct = respx{imd}(:,:,cellids(ct));
%     [N,~] = size(respct);
%     subplot(1,4,i), plot(mean(respct,2))
%     hold on
%     plot(mean(respcf,2),'r')
%     title(['GC(' num2str(cellids(cf)) '\rightarrow' num2str(cellids(ct)) ')'])
%     xlim([1,N])
%     ylim([-0.7 +0.7])
% end
% 



%% PreProcessing: Check response traces
%% Plot Hit/Miss Trial Activties of specific cell
% cc = 1;
% imd = 1;
% respcc = squeeze(resp{imd}(:,:,cc));
% 
% figure
% plot(respcc(:,:))
% hold on
% plot(mean(respcc,2),'r','linewidth',2)
% axis tight
% title(['Response during Hits for All Trials cell #' num2str(cc)])

%%% Plot Hit/Miss Activties of all cells within trial
% rr = 15;
% imd = 1;
% resprr = squeeze(respx{imd}(:,rr,:));
% %resprr = squeeze(resp{imd}(:,rr,cellids));
% 
% figure
% plot(resprr(:,:))
% axis tight
% title(['Response during Hits for All cells within trial #' num2str(rr)])

% %%
% Nc = 6;
% figure
% for i = 1:Nc
%     subplot(Nc,1,i), plot(respcc(:,i))
%     axis tight
% end
end 
