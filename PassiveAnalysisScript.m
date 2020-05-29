
   prev_young_list = {'\\Vault3\Data\Kelson\Analyzed\b233\18-07-05\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\l141\2018-05-16\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\l242\18-07-06\ToneNoiseNaive';...
    '\\Vault3\Data\Kelson\Analyzed\n226\18-07-04\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\n233\18-07-04\ToneNoise';...
    '\\Vault3\Data\Kelson\Analyzed\n242\18-07-06\ToneNoiseNaive'};


All_FRA_list =dir('//vault3/data/Kelson/analyzed/**/fra*/Fluorescence.mat')
All_FRA_list = unique({All_FRA_list.folder})';

FRA_list_noise = All_FRA_list(cellfun(@(X) ~isempty(X),...
                                                    regexp(All_FRA_list,...
                                                    '[Nn]oise')))
                                                
FRA_list_quiet =  All_FRA_list(cellfun(@(X) ~isempty(X),...
                                                    regexp(All_FRA_list,...
                                                    '[Qu]iet')))                                             

no_Fluoro = dir('//vault3/data/Kelson/analyzed/**/fra*/greenchannelregistered.raw');
no_Fluoro = {no_Fluoro.folder}';

no_Fluoro =  setdiff(no_Fluoro,All_FRA_list);

young_exp = '(\d\d\d[a-zA-Z])|([a-zA-Z]\d\d\d)';

old_exp = '[Ii][Aa]';



old_quiet = FRA_list_quiet(cellfun(@(X) ~isempty(X),...
                                                regexp(FRA_list_quiet,...
                                                old_exp)));

young_quiet = FRA_list_quiet(cellfun(@(X) ~isempty(X),...
                                                    regexp(FRA_list_quiet,...
                                                   young_exp)));


old_noise = FRA_list_noise(cellfun(@(X) ~isempty(X),...
                                                regexp(FRA_list_noise,...
                                                old_exp)));

young_noise  = FRA_list_noise(cellfun(@(X) ~isempty(X),...
                                                    regexp(FRA_list_noise,...
                                                   young_exp)));
                                               
young_noise = young_noise(cellfun(@(X) ~isempty(X),...
                                                    regexp(young_noise,...
                                                    '\\2020')));
                                                
young_noise = cat(1,prev_young_list,young_noise);       


% DataDir Name = Passive2  <- find and replace with name of DataDir
 FRA_young_noise = Fluoro_to_Table(young_noise);
 FRA_young_quiet = Fluoro_to_Table(young_quiet);
 FRA_old_noise = Fluoro_to_Table(old_noise);
 FRA_old_quiet = Fluoro_to_Table(old_quiet);
 

 
 

Passive_Paper_Analysis(FRA_young_quiet,'FRA_Quiet_Young')
Passive_Paper_Analysis(FRA_old_quiet,'FRA_Quiet_Aging')
% bayes 

