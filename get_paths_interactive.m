function DataDir = get_paths_interactive(in_path)
% finds multiple directories ineractively using UIgetfile
%

if ~exist('in_path','var')
    in_path = '\\Vault3\Data\Kelson\Analyzed'; 
end 

run = 1 
expt_id = 1;
while run
    
%%  ask user to continue
if expt_id > 1
    cnt_flg = questdlg('Continue?','Yes','No');

    if  ~ strcmp(cnt_flg,'Yes')
        run = 0; % technically not needed but for safety 
        break
    end
end


%%   
 
    
[Main_path]  = uigetdir(in_path);

if Main_path == 0 
    continue
end 

DataDir{expt_id} = Main_path
expt_id = expt_id+1; 

end 