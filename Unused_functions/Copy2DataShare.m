function Copy2DataShare(input)
%Copy2DataShare copies fluorescence traces and corresponding Psignal
%matrices to the data share
%Nikolas Francis 2016
bb=strsep(input.path,'\');
datapath = [input.savepath input.animal '\' bb{5} '\' bb{6} '\'];
sharepath = [input.sharepath input.animal '\'];
dirs = dir(datapath);
for e = 1:length(input.expname)
    for i = 1:length(dirs)
        if strcmpi(dirs(i).name,input.expname{e})
            if ~exist(sharepath)
                mkdir(sharepath);
            end
            sharedfile = [sharepath bb{5} '.' bb{6} '.' input.expname{e} '.Fluorescence.mat'];
            copyfile([datapath dirs(i).name '\Fluorescence.mat'], sharedfile)
            sharedfile = [sharepath bb{5} '.' bb{6} '.' input.expname{e} '.PsignalMatrix.mat'];
            copyfile([datapath input.expname{e} '\PsignalMatrix.mat'], sharedfile);
        end
    end
end

