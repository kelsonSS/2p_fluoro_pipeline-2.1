function Out = Fluoro_to_Table_interactive(Main_path)
% takes a cell of filepaths containing both an output file and a
% psignal_file and creates two tables from them: a table cotaining all
% neurons from that experiment sorted by freqeuency and level and a list 
% statistically significant neurons

% init
expt_id = 1 ;               % master experiment indexing variable 
run = 1;                    % boolean 
Vec_DFF_all = [];           % main data structure 
ActiveTable = table();      % activity data structure 
total_trials = 320;          % used for automated checking 
noise_settings = [-20 1 4];   % used for automated checking 

if ~exist('Main_path','var')
    Main_path = '\\vault2\Vault2Data\Kelson\Analyzed\';
end 
% main loop
while run
       
%[Cells_path,Main_path]  = uigetfile(Main_path,'Select Cell Definitions');
[Fluorescence_path,Main_path]  = uigetfile(Main_path,'Select Fluorescence');
[Psignal_file, Main_path] = uigetfile(Main_path,'Select Psignal File');
[Cell_file, Main_path] = uigetfile(Main_path,'Select Cell ID');
%[Raw_file,Main_path]  = uigetfile(Main_path,'Select registered IMG');
Paths{expt_id,1}  = Psignal_file;

    

%Cell ID handling
CellID{expt_id} = load(fullfile(Main_path,Cell_file));

% Psignal Handling
load([Main_path Fluorescence_path])
handles = WF_getPsignalInfo([Main_path Psignal_file]);


FCellCorrected = Output.FCellCorrected;
trialdur = size(FCellCorrected,1);
Neurons  = size(FCellCorrected,3);


 
 % Psignal matrix parsing
uFreqs  = handles.uFreqs ;
  uF    = length(uFreqs);
uLevels = handles.uLevels;
  uL    = length(uLevels);
  % strip non-numeric and convert to num
uLevels = cellfun(@(x) str2double(x([regexp(x,'[-0-9]')])),uLevels);
uLevels = sort(uLevels,'descend');
Freqs   = handles.Freqs  ;
Levels  = cellfun(@(x) str2double(x([regexp(x,'[-0-9]')])),handles.Levels);
FreqLevelOrder = table(Freqs,Levels);
FreqLevels = unique(FreqLevelOrder);
FreqLevels = sortrows(FreqLevels, {'Freqs','Levels'},{'ascend','descend'});



numtrials = size(FreqLevelOrder,1);


    


%% Tests 
% ensure that trials have the expected duration 
assert(trialdur == (handles.PreStimSilence + handles.PrimaryDuration+...
           handles.PostStimSilence) * handles.pfs);

     
       
  %% activity analysis
  % if we take a ttest of the mean of activity before the stimulus and
  % during the stimulus will it be significant 

  % define behaviorally relevant timepoints
  soundon = handles.PreStimSilence*handles.pfs;
  soundoff = soundon + handles.PrimaryDuration * handles.pfs;
  
  
%% calulate DF/F for each Freq/Level pair for all neurons



% baseline correction

B_Vec = repmat(mean(FCellCorrected(1:30,:,:)),[trialdur,1,1]);
Vec_DFF = (FCellCorrected -B_Vec)./B_Vec * 100;

% smoothing

gausfilt = fspecial('gaussian',[5,1],4);
Vec_DFF = imfilter(Vec_DFF,gausfilt);

% cut trials if too many trials occured
if size(Vec_DFF_all) ~= 0
    if numtrials > size(Vec_DFF_all,2)
        nt = size(Vec_DFF_all,2);
        Vec_DFF = Vec_DFF(:,1:nt,:);
        FreqLevelOrder = FreqLevelOrder(1:nt,:);
        Freqs = Freqs(1:nt);
        Levels = Levels(1:nt);
    end 
end 
%sorting 

[FreqLevelOrder, fl_ind]= sortrows(FreqLevelOrder, {'Freqs','Levels'},{'Ascend','Descend'});
 Vec_DFF = Vec_DFF(:,fl_ind,:);


  
  active = zeros(Neurons,1);
   for ii = 1:Neurons
  N = Vec_DFF(:,:,ii);
  N=N';
  
  
  [~,p]=ttest2(mean(N(:,1:soundon),2),mean(N(:,soundon+1:soundoff),2));
  active(ii,1) = p;
   end 
   % find significantly active neurons  number indicates the number of
   % stars that would be added to the figure. Zero indicates inactivity
   % while anything above 1 is active
   %
   active(active <.001) = 3;  
   active(active <.01)  = 2;  
   active(active <=.05) = 1;  
   active(active < 1  ) = 0;     
     
   %% anova based activity
   % use the fact that we presorted vec_DFF to reshape by condition
   % 10 repeats per trials 32 conditions (4 levels, 8 freqs)
%    DFF_conditions = reshape(Vec_DFF,150,10,32,[])
%    DFF_conditions_mu = squeeze(mean(DFF_conditions(soundon:soundoff,:,:,:)));
%   
%    
%    % initialize anova based signifiance list 
%    active2 = zeros(1,Neurons);
%    for ii = 1:Neurons 
%        active2(1,ii) = anova1(DFF_conditions_mu(:,:,ii)',[],'off');
%    end 
%    
%    active2_idx = active2<.01;
%    
%    DFF_cleaned = Vec_DFF(:,:,active2_idx);
   
%% cleanup 
 % create activetable 

      
  n_start = height(ActiveTable)+ 1 ;
  n_end = n_start + Neurons-1;
  % add to existing table
  temp_table = table([n_start:n_end]',active,'VariableNames',{'Neuron' 'Activity'});
 if ~exist('ActiveTable','var')
  ActiveTable = temp_table;
 else
  ActiveTable = [ActiveTable;temp_table];
  clear temp_table
 end
 
% append Vec_DFF
try
Vec_DFF_all = cat(3,Vec_DFF_all,Vec_DFF);   
catch
   warning('different filesizes detected! Skipping Experiment...')
   continue
end 


% DFF_idx = 1; 
% DFFtemp = zeros(height(FreqLevels),Neurons);
 for kk = 1:height(FreqLevels)
%           
%               disp(' ')
%              disp(strcat('Filtering: ', ...
%                   string(FreqLevels{kk,2}), ' dB; ',...
%                   string(FreqLevels{kk,1}), ' Hz'))
%               % find members of each Freq/Level to average
%               idx =ismember(Freqs, FreqLevels{kk,1}) &...
%                    ismember(Levels,FreqLevels{kk,2});
%               
%               n = sum(idx);
%               if n > 0
% %                   Ftemp = Fluoro(:,idx);
% %                   F =  nanmean(Ftemp,2)';
% %                   
% %                   B =  nanmean(F(1:handles.PreStimSilence * handles.pfs));                    
% %                   DF =( F - B)./B  ;
%                  % DF = imfilter(DF, gausfilt);
%                   DFFtemp(DFF_idx,:) = squeeze(...
%                                       nanmean(nanmean(Vec_DFF(60:100,idx,:))));
%                   
%                  
%                   
%                 
%                   DFF_idx = DFF_idx + 1 ; 
%               else
%                   continue
%               end
%               
  end
%  
%   if ~ exist('corr_idx','var')
%       corr_idx = 1
%   end 
%   
%   if ~ exist('s_corr_by_level','var')
%       s_corr_by_level = cell(4,3);
 %  end
   
%  
% s_cov = cov(DFFtemp)
% s_corr_all{corr_idx} = getCorrFromCov(s_cov) ;
% corr_idx = corr_idx +1;
% 
% HighLowIdx = true(Neurons,3);
% HighLowIdx(1:floor(Neurons*4/5),2) = false; % top 20%
% HighLowIdx(floor(Neurons/5):end,3) = false;  % bottom 20%
% for jj = 1:3
% for ii = 1:4
%     idx = ismember(FreqLevels{:,2},uLevels(ii))
%     s_cov_temp = cov(DFFtemp(idx,:));
%     s_corr_temp =  getCorrFromCov(s_cov_temp);
%     s_corr_by_level{ii,jj} = vertcat(s_corr_by_level{ii,jj},s_corr_temp(:));
% end 
% end 


clear DFFtemp
















%% ask user to continue
c_flag = questdlg('continue','yes','no');

if  ~strcmp(c_flag,'Yes') 
    run = 0;

end 
   
%append_expt_list
exptlist(n_start:n_end) = expt_id;
expt_id = expt_id + 1;


%% Get psignal and noise correlations 
 


end 
Out.DFF = Vec_DFF_all;
Out.active = ActiveTable;
Out.FreqLevelOrder = FreqLevelOrder;
Out.experiment_list = exptlist;
Out.CellID = CellID
%Out.SignalCorrByLevel = s_corr_by_level(:,1);
%Out.SignalCorrHigh = s_corr_by_level(:,2);
%Out.SignalCorrLow  = s_corr_by_level(:,3);
%Out.SignalCorr = s_corr_all;


    
