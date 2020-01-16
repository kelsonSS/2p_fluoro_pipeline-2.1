
function Out = averageBF(dataDir)

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
    
 
       
            % Passive-SNR- Tone Responsive Neurons
            DFF_mu = squeeze(nanmean(Fluoro,2));
            [~,onsets]  = max(DFF_mu);
            
            Tone_idx  = onsets > 60 & onsets <= 90;
            Noise_idx  = onsets > 30 & onsets <= 60;
            Offset_idx = onsets > 90 & onsets <= 120;
            
            idx_lst =  {Tone_idx ; Noise_idx ; Offset_idx };
            
               
           Ftemp = reshape(Ftemp, 150,8,4,10,[]);
           Ftemp = squeeze(mean(Ftemp(30:90,:,:,:,:))); 
           Ftemp = squeeze(mean(Ftemp,3));
           
           
           for ii = 1:3 
             DFF = Ftemp(:,:,idx_list{ii});
             Out{ii}.Levels = DFF;
             Out{ii}.Levels_avg = squeeze(mean(DFF,3));
           end 
             
             
             
             
               