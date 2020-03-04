
   prev_young_list = {'\\Vault3\Data\Kelson\Analyzed\b233\18-07-05\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\l141\2018-05-16\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\l242\18-07-06\ToneNoiseNaive';...
    '\\Vault3\Data\Kelson\Analyzed\n226\18-07-04\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\n233\18-07-04\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\n242\18-07-06\ToneNoiseNaive'};


All_FRA_list =dir('//vault3/data/Kelson/analyzed/**/FRA*/Fluorescence.mat');
All_FRA_list = unique({All_FRA_list.folder})';

no_Fluoro = dir('//vault3/data/Kelson/analyzed/**/FRA*/greenchannelregistered.raw');
no_Fluoro = {no_Fluoro.folder}';

no_Fluoro =  setdiff(no_Fluoro,All_FRA_list);

young_exp = 'r(\d\d\d[a-zA-Z])|([a-zA-Z]\d\d\d)';
old_exp = '[Ii][Aa]';

old_list = All_FRA_list(cellfun(@(X) ~isempty(X),...
                                                regexp(All_FRA_list,...
                                                old_exp)));

young_list  = All_FRA_list(cellfun(@(X) ~isempty(X),...
                                                    regexp(All_FRA_list,...
                                                   young_exp)));
                                               
young_list = cat(1,prev_young_list,young_list);                                               
% DataDir Name = Passive2  <- find and replace with name of DataDir



BD = BandwidthAnalysis(Passive2,'RFS',1);

BayesClassifierPassive(Passive2)



%Passive GC
GCanalTrialsModBalanced_TN(Passive2,'SNR')
GCanalTrialsModBalanced_TN(Passive2,'Tones')
GCanalTrialsModBalanced_TN(Passive2,'Noise')
GCanalTrialsModBalanced_TN(Passive2,'Offset')
GCanalTrialsModBalanced_TN(Passive2,'Off')



% bayes 

