function Out = Fluoro_to_Table(dataDir,Behavior)
% takes a cell of filepaths containing both an output file and a
% psignal_file and creates two tables from them: a table cotaining all
% neurons from that experiment sorted by freqeuency and level and a list 
% statistically significant neurons

% Behavior - If animals are actively performing we need to change the output
% structure of Vec_DFF to be Cells instead of a 3d Vector since there will
% be a variable number of trials per experiment
%
% 1 = behavior
% 0 = passive
%
%
% if dataDir is needed try get_paths_interactive


% init

if ~exist('Behavior','var')
    Behavior = 0 ;
end 

Vec_DFF_all = [];           % main data structure 
ActiveTable = table();      % activity data structure 
total_trials = 320;         % used for automated checking
    
noise_settings = [-20 1 4];   % used for automated checking 
expt_list = [];



for expt_id = 1:length(dataDir)


%[Cells_path,Main_path]  = uigetfile(Main_path,'Select Cell Definitions');
if class(dataDir) == 'char'
    Main_path = dataDir;
else
    Main_path = dataDir{expt_id};
end 
% Paths creation
Cell_file = 'CellDefinitions.mat';
Fluo_file = 'Fluorescence.mat';

dir_t = dir([Main_path] );
dir_t = {dir_t.name};
f_bl= ~ cellfun(@isempty,(regexp( dir_t  ,'_Phys_')));

if sum(f_bl) == 0
    continue
elseif sum(f_bl) == 1
Psignal_file= dir_t{f_bl};
else 
[PsigName, PsigPath] = uigetfile(Main_path);

Psignal_file =   [PsigPath,PsigName];
end 
clear dir_t
 
%Cell ID handling
CellID{expt_id} = load(fullfile(Main_path,Cell_file));

% Psignal Handling
load(fullfile(Main_path,Fluo_file))
handles = WF_getPsignalInfo(fullfile(Main_path,Psignal_file));

handles_all(expt_id) = handles; 
    

FCellCorrected = Output.FCellCorrected;
trialdur = size(FCellCorrected,1);
Neurons  = size(FCellCorrected,3);


 
 % Psignal matrix parsing
uFreqs  = handles.uFreqs ;
  uF    = length(uFreqs);
uLevels = handles.uLevels;
  uL    = length(uLevels);

  % strip non-numeric and convert to num
% uLevels = cellfun(@(x) str2double(x([regexp(x,'[-0-9]')])),uLevels);
uLevels = sort(uLevels,'descend');
Freqs   = handles.Freqs  ;
Levels  = handles.Levels;
Levels(Levels>100) = inf;
FreqLevelOrder = table(Freqs,Levels);
FreqLevels = unique(FreqLevelOrder);
FreqLevels = sortrows(FreqLevels, {'Freqs','Levels'},{'ascend','descend'});



numtrials = size(FreqLevelOrder,1);
numreps = floor(numtrials / uL / uF);

    
% red channel handling 
RedCells_file = 'RedNeuronNumber.mat'
if exist(fullfile(Main_path,RedCells_file),'file')
  RedCells = load(fullfile(Main_path,RedCells_file))
  RedCells = RedCells.n_neurons;
else 
    RedCells  = 0;
end 

if ~exist('Vec_DFF','var')
    RedCellID = [1:RedCells];
else 
     TotalNeurons = size(Vec_DFF,3) ;
    RedCellID = [RedCellID, (1:RedCells)+TotalNeurons];
end
 
 

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

B_Vec = repmat(nanmean(FCellCorrected(1:30,:,:)),[trialdur,1,1]);
Vec_DFF = (FCellCorrected -B_Vec)./B_Vec * 100;

% smoothing

gausfilt = fspecial('gaussian',[5,1],4);
Vec_DFF = imfilter(Vec_DFF,gausfilt);

% cut trials if too many trials occured
if  numtrials> total_trials
    if size(Vec_DFF_all) ~= 0
        nt = size(Vec_DFF_all,2);
        Vec_DFF = Vec_DFF(:,1:nt,:);
        FreqLevelOrder = FreqLevelOrder(1:nt,:);
        Freqs = Freqs(1:nt);
        Levels = Levels(1:nt);
    else 
        nt = total_trials;
        Vec_DFF  = Vec_DFF(:,1:total_trials,:); 
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
  
  
  [~,p]=ttest2(nanmean(N(:,1:soundon),2),nanmean(N(:,soundon+1:soundoff),2));
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
   % ex. 10 repeats per trials 32 conditions (4 levels, 8 freqs)
   
   
   counts =histcounts(log(FreqLevelOrder{:,1}),uF*uL);
   
   if mean(counts) == mode(counts)
    DFF_conditions = reshape(Vec_DFF,trialdur,mode(counts),uF *uL,Neurons);
  
   else
       DFF_conditions = [];
      rep_mode = mode(counts);
      condition = 1;
    for freq = 1:uF
        for lvl = 1:uL
           FL_idx=  FreqLevelOrder{:,1} == Freqs(freq) &... 
               FreqLevelOrder{:,2} == Levels(lvl);
           
            Vec_DFF_Temp = Vec_DFF(:,FL_idx,:);
            % in case trials have different reps
            Vec_DFF_Temp = Vec_DFF_Temp(:,1:rep_mode,:);
           plot(Vec_DFF_Temp(:,:,1))
            df_by_level_temp(freq,lvl,:) = ...
            nanmean(nanmean(Vec_DFF_Temp(soundon:soundoff,:,:),2));
            
            df_by_level_offset_temp(freq,lvl,:) =...
                nanmean(nanmean(...
                Vec_DFF_Temp(soundoff:soundoff+1*handles.pfs ,:,:),2));
            
            DFF_conditions(:,:,condition,:) = Vec_DFF_Temp;
            
             clear Vec_DFF_Temp 
             condition = condition+1;
        
        end
    end 
   end 
   if ~exist('df_by_level','var')
       df_by_level = [];
       df_by_level_offset= [];
   end 
   
   df_by_level = cat(3,df_by_level,df_by_level_temp);
   df_by_level_offset = cat(3,df_by_level_offset,df_by_level_offset_temp);
   clear df_by_level_temp
   clear df_by_level_offset_temp
   
   
DFF_conditions_mu = squeeze(mean(DFF_conditions(soundon:soundoff,:,:,:)));
DFF_conditions_mu_offset = squeeze(mean(DFF_conditions(soundoff:end,:,:,:)));
%    % initialize anova based signifiance list 
    active2 = zeros(1,Neurons);
    for ii = 1:Neurons 
        active2(1,ii) = anova1(DFF_conditions_mu(:,:,ii)',[],'off');
        active_offset(1,ii)= anova1(DFF_conditions_mu_offset(:,:,ii)',[],'off');
    end 
    
    active2_idx = active2<.01;
    offset_idx = active_offset <.01;
if ~(exist('anova_idx','var'))
    anova_idx = [];
    anova_offset_idx = [];
end
    anova_idx = cat(2,anova_idx,active2_idx) ;
    anova_offset_idx = cat(2,anova_offset_idx,offset_idx) 
%    DFF_cleaned = Vec_DFF(:,:,active2_idx);
   
%% Clean-up 

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
 
% append DFF's to output list
if Behavior == 1  % active condition

    Vec_DFF_all{expt_id} = Vec_DFF;
    
else  % passive condition

Vec_DFF_all = cat(3,Vec_DFF_all,Vec_DFF);

end 
    clear DFFtemp


% append Expt_list 
new_ids = ones( 1 + n_end - n_start  ,1)* expt_id;
 expt_list = [ expt_list ; new_ids];
 
 if class(dataDir) == 'char'
     break
 end 

end 


%% Packaging 
 
Out.DFF = Vec_DFF_all;
Out.active = ActiveTable;
Out.FreqLevelOrder = FreqLevelOrder;
Out.experiment_list = expt_list;
Out.RedCellID = RedCellID;
Out.df_by_level = df_by_level;
Out.df_by_level_offset = df_by_level_offset;
Out.CellID = CellID;
Out.DataDirs = dataDir;
Out.anova = anova_idx;
Out.Offset = anova_offset_idx;
Out.handles =handles_all;


    
