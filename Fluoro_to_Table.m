function Out = Fluoro_to_Table(dataDir,Freq_cutoff)
% takes a cell of filepaths containing both an output file and a
% psignal_file and creates two tables from them: a table cotaining all
% neurons from that experiment sorted by freqeuency and level and a list 
% statistically significant neurons

%Freq_cutoff - If trying to combine files with different numbers of
%frequencies can use Freq_cuttoff to only collect frequencies below the
%Freq_cuttoff
%

%
%
% if dataDir is needed try get_paths_interactive

% init
if ~exist('Freq_cutoff','var')
    Freq_cutoff = 65000;  % largst freq played in lab is usually 64k
end 




Vec_DFF_all = [];           % main data structure 
ActiveTable = table();      % activity data structure 
expt_list = [];
curr_expt = 1;

if ischar(dataDir)
num_expts = 1;
else 
    num_expts = length(dataDir);
end 

if isstruct(dataDir)
  for ii = 1:length(dataDir)
      dataDirtemp{ii} = fullfile(dataDir(ii).folder,dataDir(ii).name);
  end
  dataDir = dataDirtemp;
  clear dataDirtemp
    
end 

handles_all ={};
planeID = [];
for expt_id = 1:num_expts


    
%[Cells_path,Main_path]  = uigetfile(Main_path,'Select Cell Definitions');
if class(dataDir) == 'char'
    Main_path = dataDir;
elseif iscell(dataDir)
    Main_path = dataDir{expt_id};
end

if isempty(Main_path)
    continue
end

% Paths creation

Fluo_file = 'Fluorescence.mat';

dir_t = dir([Main_path] );
dir_t = {dir_t.name};
f_bl= ~cellfun(@isempty,(regexp( dir_t  ,'_Phys_')));
mat_bl= ~cellfun(@isempty,(regexp( dir_t  ,'.mat')));
f_bl = f_bl & mat_bl;
if sum(f_bl) == 0
    continue
elseif sum(f_bl) == 1
Psignal_file= dir_t{f_bl};
else 
[PsigName, PsigPath] = uigetfile(Main_path,'load Psignal file');

Psignal_file =   [PsigPath,PsigName];
end 
clear dir_t

%% XML handling 
XML = get_options_from_xml(fullfile(Main_path,'Experiment.xml'));
if ~isfield(XML, 'totalZplanes')
    XML.totalZplanes = 1;
end 

%% Psignal Handling
handles = WF_getPsignalInfo(fullfile(Main_path,Psignal_file));
total_trials =  size(handles.Freqs,1);

%% Fluorescence Handling & Cell ID handling



fluoro_files = dir(fullfile(Main_path,'Fluorescence*'));
fluoro_names = {fluoro_files.name};

celldef_files = dir(fullfile(Main_path,'CellDefinitions*'));
celldef_names = {celldef_files.name};
FCellCorrected = [];

for plane_idx = 1:length(fluoro_names)
    
    load(fullfile(Main_path,fluoro_names{plane_idx}))
    % get plane number
     plane_number = strrep(fluoro_names{plane_idx}, 'FluorescenceZ','');
     plane_number = str2num(plane_number(1));
     if isempty(plane_number)
         plane_number = 1;
     end 
     
     %Process planes 
    curr_plane = Output.FCellCorrected;
    curr_plane_ID = ones(size(curr_plane,3),1) * plane_number ;
    
    
    
    planeID = cat(1,planeID, curr_plane_ID);
    FCellCorrected = cat(3,FCellCorrected,curr_plane);
    % proccess celldef 
    CellDefinitions{expt_id}(plane_idx) = load(fullfile(...
                              Main_path,celldef_names{plane_idx}));
   
    

end 

trialdur = size(FCellCorrected,1);
trials = size(FCellCorrected,2);
Neurons  = size(FCellCorrected,3);

% if trials ~= total_trials
%     continue
% end 

fprintf('Expt %d of %d: %d frames x %d trials x %d neurons \n',...
          expt_id, num_expts, trialdur,trials, Neurons); 
      

%% Tests 
% ensure that trials have the expected duration 
%  if  trialdur ~= floor( (handles.PreStimSilence + handles.PrimaryDuration+...
%            handles.PostStimSilence) * handles.pfs/ XML.totalZplanes) ;
%        continue 
%        
%  end 

if any(size(FCellCorrected) == 0)
    continue
end


   if  trials < total_trials
       total_trials = trials;
   end 
%        

handles_all{end+1} = handles; 
          
      
      
      
 %% Psignal matrix parsing
uFreqs  = handles.uFreqs ;
  uF    = length(uFreqs);
uLevels = handles.uLevels;
  uL    = length(uLevels);
Freqs   = handles.Freqs(1:total_trials);
Levels  = handles.Levels(1:total_trials);
 
  

  % strip non-numeric and convert to num
% uLevels = cellfun(@(x) str2double(x([regexp(x,'[-0-9]')])),uLevels);
uLevels = sort(uLevels,'descend');

Levels(Levels>100) = inf;
uLevels(uLevels>100) = inf;
FreqLevelOrder = table(Freqs,Levels);


if uF > 1 
[FreqLevelOrder, fl_idx]= sortrows(FreqLevelOrder, {'Freqs','Levels'},{'Ascend','Descend'});

% if  mod(total_trials, uF * uL) ~= 0 
    try
    [FreqLevelOrder,fl_idx] = FreqCheck(FreqLevelOrder,fl_idx,total_trials);
    catch
        fprintf('Bad freqs, continuing \n')
        continue
        end 
%end 

%%change if you need to subset data
freq_idx = FreqLevelOrder{:,1} < Freq_cutoff;
if any(~freq_idx) 
   uF = uF -1 ;
   uFreqs = uFreqs(1:end-1);
    FreqLevelOrder = FreqLevelOrder(freq_idx,:);
    fl_idx = fl_idx(freq_idx);
end 

FreqLevels = unique(FreqLevelOrder);
FreqLevels = sortrows(FreqLevels, {'Freqs','Levels'},{'ascend','descend'});
else
    fl_idx = 1:trials;
end 

if ~ isfield(handles,'BackgroundNoise') || handles.BackgroundNoise(1) == -99
    noise_flag = 0;
    first_sound = handles.PreStimSilence;
else 
    first_sound  = min( handles.PreStimSilence,handles.BackgroundNoise(2) );
    noise_flag = 1;
end 
numtrials = size(FreqLevelOrder,1);
numreps = floor(numtrials / uL / uF);
first_sound_time = min(handles.PreStimSilence, handles.BackgroundNoise(2));
baseline_frames = floor(first_sound* handles.pfs / XML.totalZplanes);

% behavior parsing

%% red channel handling 
% RedCells_file = 'RedNeuronNumber.mat';
% if exist(fullfile(Main_path,RedCells_file),'file')
%   RedCells = load(fullfile(Main_path,RedCells_file));
%   RedCells = RedCells.n_neurons;
% else 
%     RedCells  = 0;
% end 
% 
% if ~exist('Vec_DFF','var')
%     RedCellID = [1:RedCells];
% else 
%      TotalNeurons = size(Vec_DFF,3) ;
%     RedCellID = [RedCellID, (1:RedCells)+TotalNeurons];
% end
%  
 


 
       
  %% activity analysis
  % if we take a ttest of the mean of activity before the stimulus and
  % during the stimulus will it be significant 

  % define behaviorally relevant timepoints
  soundon = floor(handles.PreStimSilence*handles.pfs/ XML.totalZplanes);
  soundoff = floor(soundon + handles.PrimaryDuration*handles.pfs/ XML.totalZplanes);
  
  
  
  
%% calulate DF/F for each Freq/Level pair for all neurons


% baseline correction

B_Vec = repmat(nanmean(FCellCorrected(1:baseline_frames,:,:)),[trialdur,1,1]);
Vec_DFF = (FCellCorrected -B_Vec)./B_Vec * 100;
%
%Sort DFF by Freq/Level Combination, shorten if necessary
Vec_DFF = Vec_DFF(:,fl_idx,:);

% smoothing

gausfilt = fspecial('gaussian',[5,1],4);
Vec_DFF = imfilter(Vec_DFF,gausfilt);

% cut trials if too many trials occured
% if  numtrials > total_trials
%     if size(Vec_DFF_all) ~= 0
%         nt = size(Vec_DFF_all,2);
%     else 
%         nt = total_trials;
%         Vec_DFF  = Vec_DFF(:,1:total_trials,:); 
%         FreqLevelOrder = FreqLevelOrder(1:nt,:);
%         Freqs = Freqs(1:nt);
%         Levels = Levels(1:nt);
% end 
% end  


%fl_idx2 = FreqLevelOrder.Levels < 99;
%Vec_DFF = Vec_DFF(:,fl_idx2,:);
%FreqLevelOrder = FreqLevelOrder(fl_idx2,:);


  active = zeros(Neurons,1);
  if noise_flag
      test_end = size(Vec_DFF,1);
  else 
      test_end = soundoff;
  end 
   for ii = 1:Neurons
  N = Vec_DFF(:,:,ii);
  N=N';
  
 
  [~,p]=ttest2(nanmean(N(:,1:baseline_frames),2),nanmean(N(:,soundon+1:test_end),2));
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
     
%    %% old based activity
%    % use the fact that we presorted vec_DFF to reshape by condition
%    % ex. 10 repeats per trials 32 conditions (4 levels, 8 freqs)
%    
   
  
   FLO =  table2array(FreqLevelOrder);
   FLO = sum(FLO,2);
   counts =histcounts(FLO,unique(FLO)+1); 
   
       DFF_conditions = [];
      rep_mode = mode(counts);
      if rep_mode > 10
          rep_mode = 10;
      end 
      condition = 1;
    for lvl = 1:uL
        for freq = 1:uF
         print_len = fprintf('Level %d Freq %d \n', uLevels(lvl),uFreqs(freq));
           FL_idx=  FreqLevelOrder{:,1} == uFreqs(freq) &... 
           FreqLevelOrder{:,2} == uLevels(lvl);
           
            Vec_DFF_Temp = Vec_DFF(:,FL_idx,:);
            % in case trials have different reps
            try
            Vec_DFF_Temp = Vec_DFF_Temp(:,1:rep_mode,:);
            catch
            end
          %  plot(Vec_DFF_Temp(:,:,1))
            baseline = squeeze(nanmean(Vec_DFF_Temp(1:baseline_frames,:,:),2));
            after_onset = squeeze(nanmean(Vec_DFF_Temp(soundon:soundoff,:,:),2));
            
            for nn = 1:size(Vec_DFF_Temp,3)
                df_by_level_sig_temp(lvl,freq,nn) = ttest2(baseline(:,nn),after_onset(:,nn),'alpha',.001);
            end
            
             df_by_level_temp(lvl,freq,:) = ...
             nanmean(after_onset);
         
             df_by_level_offset_temp(lvl,freq,:) =...
                  nanmean(nanmean(...
                  Vec_DFF_Temp(soundoff:soundoff+1*handles.pfs ,:,:),2));
              
             %DFF_conditions(:,:,condition,:) = Vec_DFF_Temp;  
            condition = condition+1;
            fprintf(repmat('\b',1,print_len))
         clear Vec_DFF_Temp 
        end
    end 
   
   if ~exist('df_by_level','var')
       df_by_level = [];
       df_by_level_sig = [];
       df_by_level_offset= [];
       responsive_neurons = [];
       
   end 
   
   
  
    df_by_level{expt_id} = df_by_level_temp;
    df_by_level_sig{expt_id} = df_by_level_sig_temp;
    df_by_level_offset{expt_id} = df_by_level_offset_temp;
    
    n_active_levels  =   squeeze( sum(...
                                    sum(...
                                        df_by_level_sig_temp,...
                                        1),...
                                        2) ) ; 
    responsive_neurons_expt = n_active_levels >=1;
    responsive_neurons = cat(1,responsive_neurons,responsive_neurons_expt);

    clear df_by_level_temp
    clear df_by_level_sig_temp
    clear df_by_level_offset_temp
    clear n_active_levels
    clear responsive_neurons_expt
    


%    
%    
%  DFF_conditions_mu = squeeze(mean(DFF_conditions(soundon:soundoff,:,:,:)));
%  DFF_conditions_mu_offset = squeeze(mean(DFF_conditions(soundoff:end,:,:,:)));
%  %    % initialize old based signifiance list 
%      active2 = zeros(1,Neurons);
%      active_offset = zeros(1,Neurons);
%      for ii = 1:Neurons 
%          active2(1,ii) = anova1(DFF_conditions_mu(:,:,ii)',[],'off');
%          active_offset(1,ii)= anova1(DFF_conditions_mu_offset(:,:,ii)',[],'off');
%      end 
%      
%      active2_idx = active2<.01;
%      offset_idx = active_offset <.01;
%  if ~(exist('anova_idx','var'))
%      anova_idx = [];
%      anova_offset_idx = [];
%  end
%     anova_idx = cat(2,anova_idx,active2_idx) ;
    % anova_offset_idx = cat(2,anova_offset_idx,offset_idx); 
    % DFF_cleaned = Vec_DFF(:,:,active2_idx);
  
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
 
 

if expt_id == 1 
    rep_mode_all = rep_mode;
end 

if rep_mode ~= rep_mode_all
    
    
    if rep_mode > rep_mode_all
        % too many repeats
        % get first rep_mode_all trials of each unique stimulus
        trial_idx = [1:size(Vec_DFF,2)];
        sub_idx = trial_idx(mod(trial_idx-1,rep_mode)>=rep_mode_all);
        
        Vec_DFF = Vec_DFF(:,sub_idx,:);
    elseif rep_mode < rep_mode_all
        % too few repeats
        
        if length(unique(counts)) == 1
            
            % create output with the correct number of trials
            Vec_DFF_temp = nan(size(Vec_DFF_all,1),size(Vec_DFF_all,2),size(Vec_DFF,3));
            
            
            temp_idx = 1;
            old_idx = 1;
            for ii = 1:height(unique(FreqLevels))
                
                Vec_DFF_temp(:,temp_idx:temp_idx + rep_mode - 1) = ...
                    Vec_DFF(:,old_idx:old_idx + rep_mode - 1);
                
                fprintf('old: %d-%d  new: %d-%d \n', old_idx,  old_idx + rep_mode-1,...
                    temp_idx, temp_idx + rep_mode_all-1)
                temp_idx = temp_idx + rep_mode_all;
                old_idx = old_idx + rep_mode ;
                
            end
            
            if size(Vec_DFF_temp,2) == size(Vec_DFF_all,2)
                Vec_DFF = Vec_DFF_temp;
            end
            clear Vec_DFF_temp
            
            
        end
    end
end 

% append DFF's to output list
try
Vec_DFF_all = cat(3,Vec_DFF_all,Vec_DFF);
catch
     fprintf('too few freqs, continuing \n')
     
    %continue
end 


    clear DFFtemp

    
 DFF_Z = squeeze( ( Vec_DFF_all - nanmean(nanmean(Vec_DFF_all )) )...
                  ./ nanstd(nanstd(Vec_DFF_all)));
 % baseline correction
 try
 B_DFF_Z = repmat(nanmean(DFF_Z(1:baseline_frames,:,:)),[trialdur,1,1]);  
 DFF_Z =  DFF_Z - B_DFF_Z ;       
 catch 
     continue
 end 
 
 % clean_index is defined by  the absolute variance of the first second
 % before stimulation being relatively low. this may break if prestim
 % silence is less than 1. a noisy trial often has 10-100times the variance
 % of a normal trial which makes it easy to spot in standard deviation
 % before sound onset
 Clean_idx_1 = squeeze(max(abs(nanstd(Vec_DFF_all(1:baseline_frames-1,:,:),[],2))) < 20 );
 
  Clean_idx_2 = squeeze(max(abs(nanmean(DFF_Z(1:baseline_frames-1,:,:),2))) < 2 );
 
  Clean_idx = Clean_idx_1 & Clean_idx_2;
 DFF2 = nanmean(Vec_DFF_all,2);

DFF_ab_max = max(abs(DFF2));
DFF_normalized = Vec_DFF_all./DFF_ab_max;

% append Expt_list 
new_ids = ones( 1 + n_end - n_start  ,1)* curr_expt;
assert(length(new_ids) == Neurons);
expt_list = [ expt_list ; new_ids];

curr_expt = curr_expt+1;
 
 % append datapaths 
 if ~exist('good_dirs','var')
     good_dirs = {};
 end 
  if  num_expts ~= 1
    good_dirs(1,end+1) = dataDir(expt_id);
  else 
      good_dirs{1} = dataDir;
  end 
 
% create class_Idx
 Classes = {'Inhbited','Noise','Tone_on','Tone_off','Noise_off'};
 [~,onsets]  = max(squeeze(DFF2));
% 
 Inhibited_idx = onsets <=30;
 Noise_idx = (onsets > 30 & onsets <= 60) * 2;
 Tone_idx  = (onsets > 60 & onsets <= 90) * 3;
 Offset_idx = (onsets > 90 & onsets<= 120) * 4;
 Off_idx = (onsets > 120) * 5 ;
% 
 Class_idx =  Inhibited_idx + Noise_idx + Tone_idx + Offset_idx + Off_idx;
 Class_idx = Class_idx'; 

 

 fprintf('%d neurons, %d expt-list,\n', length(expt_list),size(Vec_DFF_all,3))

 

end 


%% Packaging 

try
Out.DFF = Vec_DFF_all;
Out.DFF_Z = DFF_Z;
Out.Clean_idx= Clean_idx;
Out.DFF_norm = DFF_normalized;
Out.active = ActiveTable;
Out.responsive = responsive_neurons;
Out.FreqLevelOrder = FreqLevelOrder;
Out.experiment_list = expt_list;
Out.PlaneID = planeID;
%Out.RedCellID = RedCellID;
Out.df_by_level = df_by_level;
Out.df_by_level_sig = df_by_level_sig;
Out.df_by_level_offset = df_by_level_offset;
Out.CellID = CellDefinitions;
Out.DataDirs = good_dirs;
Out.nplanes = XML.totalZplanes;
Out.Class_idx = Class_idx;
Out.Classes = Classes;
%Out.anova = anova_idx;
%Out.Offset = anova_offset_idx;
Out.handles =handles_all;
catch
   Out = []
end 
    
function [FLO_fixed, fl_idx_fixed] = FreqCheck(FreqLevelOrder,fl_idx,total_trials)
    
   uF = unique(FreqLevelOrder{:,1});
   if length(uF) == 1 
       FLO_fixed = FreqLevelOrder;
       fl_idx_fixed = fl_idx;
       return
   end 
   
   fl_idx_fixed = [];
   FLO_fixed = table();
  
   FLO =  table2array(FreqLevelOrder);
   FLO(isinf(FLO)) = 99;
   FLO = sum(FLO,2);
   counts =histcounts(FLO,[unique(FLO);max(FLO)+1]);
   rep_mode = min(counts);
   
   remainder = 0;
   idx_pointer = 1;
   for ii = 1:length(counts)
       start_pointer = idx_pointer ;
       end_pointer = idx_pointer+rep_mode-1;
       
       fl_idx_fixed = cat(1,fl_idx_fixed,...
                           fl_idx(start_pointer:end_pointer) );
       
       FLO_fixed = cat(1,FLO_fixed,FreqLevelOrder(start_pointer:end_pointer,:));               
                       
       remainder = counts(ii) - rep_mode;
       idx_pointer = idx_pointer + rep_mode + remainder;
       
       
   end 
   
   

  
  
  
  
  





