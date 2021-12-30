function BF_Levels =  getBestFrequency(TN)

BF_Levels = squeeze(max(TN.df_by_level,[],2));


