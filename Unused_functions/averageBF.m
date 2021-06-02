
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
            
            
            
            
            
           
           for ii = 1:3 
             DFF = Fluoro(:,:,idx_lst{ii});
             %DF_by_level
             df_by_level = FRA_from_DFF(DFF,FLO);
             Out{iif,ii}.Lvls = df_by_level';
             
             % entropy 
             df_by_level=df_by_level';
             m = size(df_by_level,2);
             H = zeros(m,1);
             for nn = 1:size(df_by_level,2)
              
                H(nn) = myEntropy(df_by_level(:,nn));
             end 
             
             % packing 
             Out{iif,ii}.entropy = H;         
             Out{iif,ii}.Lvls_mu = squeeze(mean(df_by_level));
           end 
end            
             
             
             
               