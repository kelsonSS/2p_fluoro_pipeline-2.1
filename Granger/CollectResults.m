function  [results,dataDir] = CollectResults(dataDir)

% file setup
if ~exist('dataDir','var')
 gfile_c = 1 ; % gfile counter
    dataDir = {};
    dir_data = 'G:\My Drive\GrangerResults'
    while true
        dataDir{end+1,1} = uigetdir(dir_data)
        an =  input('Continue? [1/0]');
        if an == 1
            continue
        else
            break
        end
    end
end 


InsufficientCells = {};
 % file loading
for ii =  1:length(dataDir) 
    
    try
    fl = dir(fullfile(dataDir{ii},'data*')); 
    GC_Results = load(fullfile(fl(1).folder , fl(1).name));
    GC_Stats = load(fullfile(fl(2).folder , fl(2).name));
    catch 
    InsufficientCells{end+1} = dataDir{ii}; 
        continue
    end 
    
     % GC strength is the mean of the absolute value of all nonzero J
     % statistics 
     
    GCJL = num2cell(GC_Results.GCJL,[1,2]);
    GCJL = squeeze(...
               cellfun(@(x) mean(abs( x(x ~= 0))),GCJL,'UniformOutput',1));
     
    GCJH = num2cell(GC_Results.GCJH,[1,2]);
    GCJH = squeeze(...
               cellfun(@(x) mean(abs( x(x ~= 0))),GCJH,'UniformOutput',1));    
    
    
    if ~exist('results','var')
        results = GC_Stats;
        results.GCstrength_L = GCJL;
        results.GCstrength_H = GCJH;
    else
        results.GCnumbers(:,end+1) = GC_Stats.GCnumbers;
        results.GClengths(:,end+1) = GC_Stats.GClengths;
        results.GCangles(:,end+1) = GC_Stats.GCangles;
        results.GCavgL(:,end+1) = GC_Stats.GCavgL;
        results.GCstrength_H(:,end+1) = GCJH;
        results.GCstrength_L(:,end+1) = GCJL;
           
    end  
end

% correcting orientation

 results.Exptlist = dataDir;
 results.GCnumbers = results.GCnumbers';
 results.GClengths = results.GClengths';
 results.GCangles = results.GCangles';
 results.GCavgL = results.GCavgL';
 results.GCstrength_H = results.GCstrength_H';
 results.GCstrength_L = results.GCstrength_L';

 
 %% analysis  and cleanup
results.GCnumbers_mu = mean(results.GCnumbers);
results.GCnumbers_std = std(results.GCnumbers);
results.InsufficientCells = InsufficientCells;


