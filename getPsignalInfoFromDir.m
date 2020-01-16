function handles = getPsignalInfoFromDir(Main_path)

% this function grabs the experiment file from the directory supplied and
% then calls WF_getPsignalinfo to extrct the relevant timing parameters.

% WARNING: this function is predicated on the fact that there should only
% be one psignalfile per experiment - I.E only one psignalfile in the
% directory supplied. If there are more than one the function defaults to
% using uigetfile for the user to pick the correct file 

% - Kelson Shilling-Scrivo 2019 

dir_t = dir([Main_path] );
dir_t = {dir_t.name};
Psignal_files= ~ cellfun(@isempty,(regexp( dir_t  ,'_Phys_')));


if sum(Psignal_files) == 0
    return
elseif sum(Psignal_files) == 1
Psignal_file= dir_t{Psignal_files};
else 
[PsigName, PsigPath] = uigetfile(Main_path);

Psignal_file =   [PsigPath,PsigName];
end 


handles = WF_getPsignalInfo(fullfile(Main_path,Psignal_file));



end 