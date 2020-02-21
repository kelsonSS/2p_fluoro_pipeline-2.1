function Out = Fluoro_to_Table_interactive(in_path)
% takes a cell of filepaths containing both an output file and a
% psignal_file and creates two tables from them: a table cotaining all
% neurons from that experiment sorted by freqeuency and level and a list 
% statistically significant neurons

% Warning! - it is currently up to the user to select all signal files from
% the same experiment type


% init
expt_id = 1 ;               % master experiment indexing variable 
expt_id_flag = 0;           % flag for properly indexing
run = 1;                    % boolean 
Vec_DFF_all = [];           % main data structure 
ActiveTable = table();      % activity data structure 
total_trials = 320;          % used for automated checking 
noise_settings = [-20 1 4];   % used for automated checking 
exp_files = {};
exp_bool = [];
Vec_DFF_all= [];
dataDir = [];

if ~exist('in_path','var')
    in_path = '\\Vault3\Data\Kelson\Analyzed'; 
end 


while run

%%  ask user to continue
if expt_id > 1
    cnt_flg = questdlg('Continue?','Yes','No');

    if  ~ strcmp(cnt_flg,'Yes')
        run = 0; % technically not needed
        break
    end
end

if expt_id_flag
    expt_id = expt_id + 1 ;

    expt_id_flag = 0

end


%%   
 
    
[Main_path]  = uigetdir(in_path);

if Main_path == 0 
    continue
end 

dataDir{expt_id} = Main_path;
exp_files{expt_id} = Main_path;

% Paths creation
Cell_file = 'CellDefinitions.mat';
Fluo_file = 'Fluorescence.mat';


dir_t = dir([Main_path] );
dir_t = {dir_t.name};
f_bl= ~ cellfun(@isempty,(regexp( dir_t  ,'_Phys_')));
if sum(f_bl) == 1
Psignal_file= dir_t{f_bl};
elseif sum(f_bl) > 1
psignal_file  = uigetfile(Main_path)
else 
  continue 
end 
clear dir_t
    
%Cell ID handling
CellID{expt_id} = load(fullfile(Main_path,Cell_file));

% Psignal Handling 
 
load(fullfile(Main_path,Fluo_file))
handles = WF_getPsignalInfo(fullfile(Main_path,Psignal_file));


FCellCorrected = Output.FCellCorrected;
trialdur = size(FCellCorrected,1);
Neurons  = size(FCellCorrected,3);
exp_bool = [exp_bool, repmat(expt_id,1,Neurons)]; 

 
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
 
 

 % Psignal matrix parsing
uFreqs  = handles.uFreqs ;
  uF    = length(uFreqs);
uLevels = handles.uLevels;
  uL    = length(uLevels);
 Freqs  = handles.Freqs;
 Levels = handles.Levels;
 
  % strip non-numeric and convert to num
 if iscell(uLevels) 
    uLevels = cellfun(@(x) str2double(x([regexp(x,'[-0-9]')])),uLevels);
 end 

 if iscell(Levels) 
   Levels  = cellfun(@(x) str2double(x([regexp(x,'[-0-9]')])),handles.Levels);
 end 
 
 uLevels = sort(uLevels,'descend');
 

FreqLevelOrder = table(Freqs,Levels);
FreqLevels = unique(FreqLevelOrder);
FreqLevels = sortrows(FreqLevels, {'Freqs','Levels'},{'ascend','descend'});
numtrials = size(FreqLevelOrder,1);


    


%% Tests 
% ensure that trials have the expected duration and number of trials

 trialdur_test = trialdur == (handles.PreStimSilence + handles.PrimaryDuration+...
           handles.PostStimSilence) * handles.pfs;

if exist('Vec_DFF','var')       
    Nrep_test = numtrials >=  size(Vec_DFF,2);
else 
    Nrep_test = 1;
end 
     
if ~(trialdur_test && Nrep_test)
    continue
end 

       
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
if  any(size(Vec_DFF_all) ~= 0)
    nt = size(Vec_DFF_all,2);
    if numtrials > nt
        Vec_DFF = Vec_DFF(:,1:nt,:);
        FreqLevelOrder = FreqLevelOrder(1:nt,:);
        Freqs = Freqs(1:nt);
        Levels = Levels(1:nt);
    end
    
  if  numtrials < nt
      small_trial_text = sprintf('expected trials = %d actual = %d. add?',...
          nt, numtrials);
     small_flag = questdialog(small_trial_text, 'yes' ,'no')
       
      if strcmp(small_flag,'yes')
        trials_to_add = nt - numtrials;
        blank_data = nan( size(Vec_DFF_all,1),...
                              trials_to_add,...
                              size(Vec_DFF_all,3));
          Vec_DFF = cat(2,Vec_DFF,blank_data);                
      else 
          continue
      end 
  end 
      
end



% add trial to total 
Vec_DFF_all = cat(3,Vec_DFF_all,Vec_DFF);

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
 

clear DFFtemp



expt_id_flag = 1; 

end 
%% Packaging 
 
Out.DFF = Vec_DFF_all;
Out.Active = ActiveTable;
Out.FreqLevelOrder = FreqLevelOrder;
Out.Experiment_files = exp_files;
Out.Experiment_list = exp_bool;
Out.CellID = CellID;
Out.RedCells = RedCellID;
Out.DataDirs = dataDir;
Out.Handles =  handles;   




    
