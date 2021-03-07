function needs_registration =CreateRegistrationList(savepath)

if ~exist('savepath','var')
    savepath = '\\vault3\data\kelson\analyzed';
end 


 extracted = dir( fullfile(savepath, '\**\GreenChannel*.raw'));
 extracted = struct2cell(extracted);
 extracted = extracted(2,:)';
 
 registered = dir( [savepath, '\**\GreenChannelRegistered*.raw']);
 registered = struct2cell(registered);
 registered = registered(2,:)';
 
 
 needs_registration =  setdiff(extracted,registered);
