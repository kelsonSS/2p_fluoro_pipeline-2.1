%% Cell_ID script
% using Passive.DataDirs which is a list of folders containing all experimental data
% This script will extract which cells were used for GC_analysis and then 
% uses that to plot the FRAs of those neurons
% 
Outpath = 'C:\Users\Kelson\Dropbox\Passive_paper\';

Cell_ID.All = GC_get_Cell_IDs(Passive.DataDirs,'SNR'); % ALL Data 
Cell_ID.Tones = GC_get_Cell_IDs(Passive.DataDirs,'Tones');
Cell_ID.Noise = GC_get_Cell_IDs(Passive.DataDirs,'Noise');
Cell_ID.Offset = GC_get_Cell_IDs(Passive.DataDirs,'Offset');

Cell_fields =  fieldnames(Cell_ID);

FRA =struct();
for fld = 1:length(Cell_fields)
    field = Cell_fields{fld}; 
    FRA.(field) =  FRA_from_Cell_ID(Passive.DataDirs ,Cell_ID.(field),Outpath,field);
end 