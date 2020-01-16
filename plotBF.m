function   [MuTable] = plotBF()

extra_plots = 1; % adds additional debugging/misc plots
pause_flag = 0; % pause after creation of every plot
% loading experiement info 
if ~exist('Main_path', 'var')
Main_path = '\\vault2\Vault2Data\Kelson\Analyzed\';
end 

% if ~exist('save_path',var)
%    save_path = []  
% end 

[Cells_path,Main_path]  = uigetfile(Main_path,'Select Cell Definitions'); 
[Fluorescence_path,Main_path]  = uigetfile(Main_path,'Select Fluorescence');
[Psignal_file, Main_path] = uigetfile(Main_path,'Select Psignal File');
[Raw_file,Main_path]  = uigetfile(Main_path,'Select registered IMG');


Cells = load([Main_path Cells_path]);
load([Main_path Fluorescence_path])
handles = WF_getPsignalInfo([Main_path Psignal_file]);
try
    muImg = struct2array(load([Main_path Raw_file]))';
catch
    warning('no AvgImg Found! Press any key to continue '  )
    pause
    muImg = nan(512);
end 

FCellCorrected = Output.FCellCorrected;
%clear Output 
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
trialdur2 = (handles.PreStimSilence + handles.PrimaryDuration+...
           handles.PostStimSilence) * handles.pfs;

%% Tests 
% ensure that trials have the expected duration 
try
assert(trialdur == trialdur2)
catch
   assert( ~(abs(trialdur - trialdur2) - 1))  % test if there is a 1 frame difference 
end 
       



% reps    = length(uLevels) * length(uFreqs)

% Levels(Levels == uLevels(1)) = 65;
% Levels(Levels == uLevels(2)) = 55;
% Levels(Levels == uLevels(3)) = 45;


% FreqLevels(:,1) = repmat(unique(Freqs),size( unique(Levels) ));
% FreqLevels(:,2) = nan;
% 
%  m = length(unique(Freqs));
%  mapping = uLevels;
%  id = 1;
% 
%
%  
%   for ii = 1:length(mapping)  
%      
%      FreqLevels( id: id + m-1 ,2 )= mapping( ii );
%      id= id + m;
%   end

  %% activity analysis
  % if we take a ttest of the mean of activity before the stimulus and
  % during the stimulus will it be significant 

  % define behaviorally relevant timepoints
  soundon = handles.PreStimSilence*handles.pfs;
  soundoff = soundon + handles.PrimaryDuration * handles.pfs;
  
  
  
  
   for ii = 1:Neurons
  N = FCellCorrected(:,:,ii);
  N=N';
  
  Silence = mean(N(:,1:Prestimsilence));
  
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
     
 RAN_idx =  1;  
ActiveTable = table([1:Neurons]',active,'VariableNames',{'Neuron' 'Activity'});
% select all active neurons and output to an array
ANeurons    = table2array(ActiveTable(ActiveTable.Activity  > 0, 1));
RANeurons   = table2array(ActiveTable(ActiveTable.Activity >= RAN_idx, 1)) ;
UNeurons    = table2array(ActiveTable(ActiveTable.Activity  == 0, 1));
 

  
  
%% calulate DF/F for each Freq/Level pair for all neurons
gausfilt = fspecial('gaussian',[20,1],5);


% baseline correction 
B_Vec = repmat(mean(FCellCorrected(1:30,:,:)),[trialdur,1,1]);
Vec_DFF = (FCellCorrected -B_Vec)./B_Vec;
% smoothing

Vec_DFF = imfilter(Vec_DFF,repmat(gausfilt),[1,trials,Neurons]);



%% Optional Conversion to Table format.  
% vec_order = repmat(FreqLevelOrder, [Neurons,1]);
% vec_neurons =repmat([1:Neurons], [numtrials,1]);
% vec_neurons = vec_neurons(:);
% 
% Ftable = [ table(vec_neurons), vec_order , table(reshape(Vec_DFF,150,[])') ];
%     
% Ftable.Properties.VariableNames = {'Neuron','Freqs','Level','DeltaF'};



%% Create a table of the mean of each Freq/Level Pair
MuTable = [];

 n_frames = size(FCellCorrected,1);
 n_trials = size(FCellCorrected,2);
 n_neurons = size(FCellCorrected,3);
 Fmu = squeeze( mean(FCellCorrected,2));


  B = mean(Fmu(1:30,:));
  B_mat = zeros(size(FCellCorrected));
  for ii = 1:n_neurons
      B_mat(:,:,ii) = B(ii);
  end 
 
  DFF_all = (FCellCorrected - B_mat)./B_mat *100;
  DFF_all = imfilter(DFF_all,gausfilt);
  clear B_mat

reps = height(FreqLevels);
for jj = 1:length(RANeurons)
      neuron_ID = RANeurons(jj);
      Neuron = repmat(neuron_ID,reps,1);
      Fluoro = FCellCorrected(:,:,neuron_ID);
     
      DFFtemp = zeros(height(FreqLevels),size(Fluoro,1));
      DFF_idx = 1 ;
      for kk = 1:height(FreqLevels)
          
              disp(' ')
             disp(strcat('Filtering Neuron: ', num2str(neuron_ID) ,';', ...
                  string(FreqLevels{kk,2}), ' dB; ',...
                  string(FreqLevels{kk,1}), ' Hz'))
              % find members of each Freq/Level to average
              idx =ismember(Freqs, FreqLevels{kk,1}) &...
                   ismember(Levels,FreqLevels{kk,2});
              
              n = sum(idx);
              if n > 0
%                   Ftemp = Fluoro(:,idx);
%                   F =  nanmean(Ftemp,2)';
%                   
%                   B =  nanmean(F(1:handles.PreStimSilence * handles.pfs));                    
%                   DF =( F - B)./B  ;
                 % DF = imfilter(DF, gausfilt);
                  DFFtemp(DFF_idx,:) = squeeze(...
                                      nanmean(DFF_all(:,idx,neuron_ID),2));
                  
                 
                  
                
                  DFF_idx = DFF_idx + 1 ; 
              else
                  continue
              end
              
      end
     
      
      t_temp = [table(Neuron), FreqLevels , table(DFFtemp)];
    
    
      MuTable = [MuTable;t_temp];
      
      
 %% Debugging Plot checking that DFF and DFF Vec give similar results
%       figure(2)
%       subplot(2,1,1)
%       plot(Vec_DFF(:,:,ii))
%       hold on 
%       plot(mean(Vec_DFF(:,:,ii),2),'k','LineWidth',2)
%       title(sprintf('neuron %d dF/F all trials',ii))
%       hold off
%        subplot(2,1,2)     
%        plot(DFF')
%        hold on 
%        plot(mean(DFF)','k','LineWidth',2);
%        title('dF/F mu by freq*level')
%        hold off
%        pause
 %%      
     
end
  
  


  MuTable.Properties.VariableNames = {'Neuron' 'Freq' 'Level' 'DeltaF'};
%   MuTable.Level = cellfun(@str2num,MuTable.Level) ;
 MuTable = sortrows(MuTable,{'Neuron','Freq','Level'},{'ascend','ascend','descend'});
 
 Atable = grpstats(MuTable(:,[1,2,4]),{'Neuron','Freq'}); 

 
 
 
 
 
 %% Plotting Mean dF/F and FRA for each active neuron 
 lvls = length(uLevels);
 frqs = length(uFreqs);
 RANs = length(RANeurons);
 
 gauss = fspecial('gaussian',[50,1]);
 df_by_level_all = zeros(lvls,frqs,RANs);
 
 
 for ii = 1:RANs
      
    % extract a single neurons Responses
    nn  = RANeurons(ii);
    df_by_level = (MuTable(MuTable.Neuron == nn,4));
    df_by_level = table2array(( df_by_level));
    df_by_level = df_by_level - nanmean(df_by_level(:,1:30),2); % set baseline to zero
    df_by_level_mu = max(df_by_level(:,soundon:soundoff),[],2);
    df_by_level_flat(:,ii)  = df_by_level_mu;
    df_by_level_all(:,:,ii) = reshape(df_by_level_mu,lvls,frqs,[]);
   
  % plot the responses 
    figure
    set(gcf,'Position',[5   129   648   837])
    subplot(1,2,1)
    mu = mean(df_by_level);
    hold on 
    plot(mu,'b')
    %plot(mean(Vec_DFF(:,:,nn),2),'k','LineWidth',2)
    aa = axis;
    plot([aa(1), aa(2)], [0 0 ] ,'--k')
    plot([30 30], [aa(3),aa(4)],'--r')
    plot([60 60], [aa(3),aa(4)],'--g')
    plot([90 90], [aa(3),aa(4)],'--g')
    plot([120 120], [aa(3),aa(4)],'--r')
    
%     subplot(2,2,2)
%     
%     myFRA(uFreqs,uLevels(1),df_by_level_all(1,:,ii),nn)
    subplot(1,2,2)
    myFRA(uFreqs,uLevels,df_by_level_all(:,:,ii),nn);
    % yticklabels({'inf','20' ,'+10','00'})
    
    
    if extra_plots == 1
        %% Mean by Freq/Level
        %  This plot shows the average response by frequency and level combination
        %  Requires df_by_level and a list of unique freqs and levels 
        
        
        %  reshape df_by_level to be sorted row major instead of colum
        %  major since that is the behavior of subplot
        
        df_by_level_xy = reshape(df_by_level,lvls,frqs,[]);
        df_by_level_xy = permute(df_by_level_xy,[2,1,3]);
        df_by_level_xy = reshape(df_by_level_xy,frqs*lvls,[]);
        plt = 1 ;
         figure
        set(gcf,'Position', [756 122 648 837])
     for lvl = 1:lvls
         for frq = 1:frqs
          % plot avg 
          subplot(lvls,frqs,plt)
          hold on
          plot(df_by_level_xy(plt,:)','r','LineWidth',2)
          axis tight 
          df_by_level_mu2(plt) = max(df_by_level_xy(plt,soundon:soundoff));
          % add point for mean 
          x_pt = floor((soundon+soundoff)/2);
         % scatter(x_pt,df_by_level_mu2(plt));
          % plot axes 
              aa = axis;
        plot([aa(1), aa(2)], [0 0 ] ,'--k')
       % plot([30 30], [aa(3),aa(4)],'--r')
       % plot([60 60], [aa(3),aa(4)],'--g')
       % plot([90 90], [aa(3),aa(4)],'--g')
       % plot([120 120], [aa(3),aa(4)],'--r')
          % titles on outside rows/columns
          if lvl == 1
          title( uFreqs(frq))
          end 
          if frq == 1
              ylabel(['dF/F / ' num2str(uLevels(lvl))])
          end
          % increment plot idx
         plt= plt+1;  
         end
         
     end 
    linkYlim(gcf)
    
    
   try
    assert(all(df_by_level_mu == df_by_level_mu2'));
    % if assertion fails uncomment this code to evaluate why 
%     figure 
%     subplot(1,2,1)
%     myFRA(uFreqs,uLevels,df_by_level_mu,nn)
%     subplot(1,2,2)
%     myFRA(uFreqs,uLevels,df_by_level_mu2',nn)
    catch
    end 
    
    end 
    
   if pause_flag
    pause
   end 
    
 end 
     
 %% Spatial Plotting of BFs
 

 % creating colormap and appending Nans for undefined neurons 
 level_colors  = [jet(uL); [nan nan nan] ];
 Freq_colors = [jet(uF) ; [nan nan nan]  ] ;
 
 sFreq_idx = zeros(RANs,1);
 sLevel_idx =zeros(RANs,1);
 
 for ii = 1:size(df_by_level_all,3)
     
     s_df_temp = df_by_level_all(:,:,ii);
     max_val  = max(s_df_temp(:));
     if max_val ~= 0 
        [s_idx] =  find(s_df_temp == max(s_df_temp(:)));
        [sLevel,sFreq] = ind2sub(size(s_df_temp),s_idx);
     
     else % non-active neurons
        sFreq = uF+1;
        sLevel= uL+1;
     end 
     
     sFreq_idx(ii) = sFreq;
     sLevel_idx(ii) = sLevel;
     
 
 
 end 
 clear df_temp 
 %FreqMap
 figure
 subplot(1,2,1)
hold on 
colormap('gray')
imagesc(muImg)
 gscatter(Cells.ptsIdx(RANeurons,3),Cells.ptsIdx(RANeurons,2),...
    uFreqs(sFreq_idx))
 title('Frequency Response')
 setScatterSettings(gcf)
%levelmap
 subplot(1,2,2) 
 hold on 
imagesc(muImg)
 gscatter(Cells.ptsIdx(RANeurons,3),Cells.ptsIdx(RANeurons,2),...
    uLevels(sLevel_idx))

title('Level Response')

setScatterSettings(gcf)
 

 
%% KMeans Clustering Noise Vs Tone Responsive neurons 
try
Vec_DFF_mu = squeeze(mean(Vec_DFF(:,:,RANeurons),2));
Vec_DFF_mu(Vec_DFF_mu<0.01) = 0;
Vec_DFF_mu = Vec_DFF_mu./max(Vec_DFF_mu);

k = 4;
Kmeans_onset = mean(Vec_DFF_mu(31:61,:));
Kmeans_tone = mean(Vec_DFF_mu(61:90,:)) ;
Kmeans_offset = mean(Vec_DFF_mu(91:120,:));
Kmeans_colors= [jet(k); [0 0 0] ]; % adding zeros to deal with Nans


Kmeans_idx =  kmeans([Kmeans_onset;Kmeans_tone]',k) ;
Kmeans_idx(isnan(Kmeans_idx)) = k+1; %Nans will now have an associated color 
figure
scatter(Kmeans_onset,Kmeans_tone,30, ... 
        Kmeans_colors(Kmeans_idx,:),'filled')
hold on 
% plot 45 degree line       
xy = linspace(0,1,20);
plot(xy,xy, 'k--') 
scatter(Kmeans_onset,Kmeans_tone)

hold off
xlabel('Normalized Mean Onset DF/F')
ylabel('Normalized Mean Tone DF/F')

% Spatially plot the BFs in space (requires Cell IDs)
figure
 scatter(Cells.ptsIdx(RANeurons,3),Cells.ptsIdx(RANeurons,2),50,...
         Kmeans_colors(Kmeans_idx,:),'filled',...
         'MarkerEdgeColor','b')
 
 hold on 
 
  % find centroids 
  Kmeans_centroid = zeros(k,2);
  for ii = 1:k
      xytemp = Cells.ptsIdx(Kmeans_idx == ii,[2 3]);
      Kmeans_centroid(ii,:) = mean(xytemp);
      scatter(Kmeans_centroid(ii,2),Kmeans_centroid(ii,1),50,'*',...
              'MarkerFaceColor', Kmeans_colors(ii,:),...
              'MarkerEdgeColor','k' ) 
          
  end 
%  
 
 clear xytemp
 hold off 

catch
end 
 %% Signal And Noise Correlations: Method 1-FRA Based 
 % 
 % 
 % plot the signal and noise correlations by first creating the covariance
 % matrix where each entry ij corresponds to the cov(I,J) and the
 % main diagonal is equal to var(I)
 % signal correlations equation = Cov(I,J)/sqrt(var(I)*var(j)) 
 %
 % s = signal n = noise 
s_cov_1 = cov(df_by_level_flat);
s_corr_1 = getCorrFromCov(s_cov_1);



% subtract the mean off each tone/level combination from each presentation


DFF_all_mu = squeeze(mean(DFF_all(:,:,RANeurons),2));
DFF_all_mu_subtracted  = nan(n_trials,RANs);

DFF_trial = squeeze(mean(DFF_all(soundon:soundoff,:,:)));
for jj = 1:length(RANeurons)
    neuron_ID = RANeurons(jj);
    
    Fluoro = DFF_all(:,:,neuron_ID);
    
    DFFtemp = zeros(height(FreqLevels),size(Fluoro,1));
    DFF_idx = 1 ;
    for kk = 1:height(FreqLevels)
        
        disp(' ')
        disp(strcat('Filtering Neuron: ', num2str(neuron_ID) ,';', ...
            string(FreqLevels{kk,2}), ' dB; ',...
            string(FreqLevels{kk,1}), ' Hz'))
        % find members of each Freq/Level to average
        idx =ismember(Freqs, FreqLevels{kk,1}) &...
            ismember(Levels,FreqLevels{kk,2});
        
        
        
        DFF_all_mu_subtracted(idx,jj) =  DFF_trial(idx,jj) -...
            df_by_level_flat(kk,jj);
        
        
    end
end 
clear DFFtemp

n_cov_1 = cov(DFF_all_mu_subtracted);
n_corr_1 = getCorrFromCov(n_cov_1);

% repeat analysis, spliting data into the top and bottom 50% of neurons by
% average DF/F

DFF_avg =  mean(df_by_level_flat);
DFF_avg(DFF_avg<0) = 0;
a= 1; % if a == 1, DFF_idx splits at mean DFF
DFF_idx = DFF_avg > a * squeeze(mean(DFF_avg)) ;


s_corr_high = getCorrFromCov(cov(df_by_level_flat(:,DFF_idx)));
s_corr_low   = getCorrFromCov(cov(df_by_level_flat(:,~DFF_idx)));

n_corr_high = getCorrFromCov(cov(DFF_all_mu_subtracted(:,DFF_idx)));
n_corr_low  = getCorrFromCov(cov(DFF_all_mu_subtracted(:,~DFF_idx)));

figure 
subplot(1,2,1);title('Signal Correlation ')
hold on 
histogram(s_corr_1)
histogram(s_corr_high)
histogram(s_corr_low)
legend('all','high','low')
subplot(1,2,2);title('Noise Correlation ')
hold on 
histogram(n_corr_1)
histogram(n_corr_high)
histogram(n_corr_low)
legend({'all','high','low'})



%% 
% PCA-based population vector
% commented out because of low interpretability. Feel Free to try out on 
% vec_DFF_RAN = Vec_DFF(:,:,RANeurons);
% % normalize responses by the max response of each neuron 
% vec_max = squeeze(max(max(vec_DFF_RAN)));
% vec_max = permute(repmat(vec_max,1,320,150),[3,2,1]);
% 
% vec_DFF_RAN = vec_DFF_RAN ./vec_max;
% 
% clear vec_max
% 
% 
% [p.coef,p.score,p.latent,p.t2,p.explained,p.mu]=...
%         pca(squeeze(mean(Vec_DFF_RAN,2)));
% 
%  %  percent explained by top eginvectors 
%     p.top2eig = sum(p.explained(1:2));
%     p.top3eig = sum(p.explained(1:3));
%  %  get top 3 axes for plotting  
%     p.axis1 = p.coef(1,:)';
%     p.axis2 = p.coef(2,:)';
%     p.axis3 = p.coef(3,:)';
%     
%     % vary the size parameters of the scatterplot
%     L=size(vec_DFF_RAN,1);
%     stim_mu = floor(mean(soundon,soundoff));
%     sz = fspecial('gaussian',[L,1],20)*1000;
%     sz = circshift(sz, stim_mu - floor(L/2) );
%     % plotting
%     Colors = parula(frqs);
%     figure
%     dx_mu = zeros(frqs,1);
%     dy_mu = zeros(frqs,1); 
%     dz_mu = zeros(frqs,1);
%     x_mu = zeros(frqs,1);
%     y_mu =zeros(frqs,1);
%     z_mu = zeros(frqs,1);
%     for trial = 1:size(Vec_DFF_RANs,2)
%         x= squeeze(vec_DFF_RAN(:,trial,:)) * p.axis1;   
%         y= squeeze(vec_DFF_RAN(:,trial,:)) * p.axis2; 
%         z= squeeze(vec_DFF_RAN(:,trial,:)) * p.axis3;
%         dx = x(2:end)-x(1:end-1);
%         dy = y(2:end)-y(1:end-1);
%         dz = z(2:end)-z(1:end-1);
%         freqColor_idx =  (uFreqs == Freqs(trial));
%         x_mu(freqColor_idx) = x_mu(freqColor_idx)+ x(1)/10;
%         y_mu(freqColor_idx) =  y_mu(freqColor_idx)+y(1)/10;
%         z_mu(freqColor_idx) =  z_mu(freqColor_idx)+z(1)/10;
%         dx_mu(freqColor_idx) = dx_mu(freqColor_idx)+dx/10;
%         dy_mu(freqColor_idx) = dy_mu(freqColor_idx)+dy/10;
%         dz_mu(freqColor_idx) = dz_mu(freqColor_idx)+dz/10;
%         subplot(1,2,1)
%         hold on 
%         quiver(0,0,dx,dy,2,'Color',Colors(freqColor_idx,:))
%         subplot(1,2,2)
%         hold on
%         quiver3(x(1:end-1),y(1:end-1),z(1:end-1),dx,dy,dz,2,'Color',Colors(freqColor_idx,:))
%             pause 
%     end 
%     figure
%     subplot(1,2,1)
%     quiver(x_mu,y_mu,dx_mu,dy_mu)
%     subplot(1,2,2)
%     quiver3(x_mu,y_mu,z_mu,dx_mu,dy_mu,dz_mu)

end 

function setScatterSettings(gcf)

ax = findall(gcf,'Type','Axes');
lgnd = findall(gcf,'Type','Legend');

set(ax,'visible','off')
set(findall(gcf,'Type','text'),'visible','on')
set(ax,'color','none')
set(lgnd,'Location','east')
lgnd_pos = get(lgnd,'position');
% shift right
try 
lgnd_new = lgnd_pos + [.1 , 0 ,0 , 0];
set(lgnd,'position',lgnd_new)
catch 
lgnd_new = cellfun(@(x) x+[.1 0 0 0] ,lgnd_pos,'UniformOutput',0);
for ii = 1:length(lgnd_pos)
    set(lgnd(ii),'position',lgnd_new{ii})
end 

end 


end 
  