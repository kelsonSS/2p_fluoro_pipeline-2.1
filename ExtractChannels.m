function ExtractChannels(input)
%ExtractChannels is an efficient vectorized program that loads the raw image
%sequence from ThorImage and extracts just the red and green channels.
%Note: does not extract more than 2 channels and assumes that the green
%channel comes before the red channel
%Note2: this method is extremely RAM intensive reccommend 24G of ram for
%10G of video.This estimate may not be accurate YMMV.If this method breaks
%try adapting ExtractRedGreenChannels for your use. 
%Kelson Shilling-Scrivo 2018
tic();
%Load image sequences into memory
bb=strsplit(input.path,'\'); %changed from strsep to strsplit
filePreparation(input)
%% Extract red and green channels, and store them in the local path
for i = 1:length(input.expname)
    disp(['Extracting ' input.path '\' input.expname{i}])
    
    newpath = fullfile(input.savepath, bb{end}, input.expname{i});
    xmlpath = fullfile(newpath,'Experiment.xml');
    opts = get_options_from_xml(xmlpath);
    
    
   
    out_path = fullfile(newpath,'greenchannel.raw');
    
    % chcek to see if the file has already be extracted 
    if ~ isempty(dir(out_path))
        newfile = dir(out_path);
        oldfile = dir(fullfile(input.path,input.expname{i},'\Image_0001_0001.raw'));         
        if newfile.bytes == oldfile.bytes 
            continue
        end
    end 
    
    if opts.numchannels == 1
        
        old_path = fullfile(input.path,input.expname{i},'Image_0001_0001.raw');
        try
        copyfile(old_path,out_path, 'f')
        catch
            continue
        end
        
    elseif opts.numchannels == 2
        
        out_path = fullfile(newpath,'greenchannel');
        img = fopen(fullfile(input.path,input.expname{i},'Image_0001_0001.raw'));
       
        try
        Green = fread(img,'uint16=>uint16');
        catch 
            fclose(img);
            continue
        end
        fclose(img);
       
        
        green_channel = opts.greenChannel   ; 
        red_channel = opts.redChannel;
        Green = reshape(Green, opts.dimX,opts.dimY,[]);
        green_L = size(Green,3);
        
        version = opts.version(1:3);
        % determine the number of channels saved 
        if version == '2.1'
            nchannels = 4 ;
            
             Red = Green(:,:,red_channel:nchannels:end);
        % remove all other channels to make green
        not_green_idx = 1:size(Green,3); 
        not_green_idx(green_channel:nchannels:end) = [];
        Green(:,:,not_green_idx) =[];
        
            
            
            
        elseif strmatch('3.1',version) || strmatch(version,'3.0')

            nchannels = 2 ;
            
          
            Red = Green(:,:,[2:2:end]);
            Green(:,:,[2:2:end]) = [];
        
        end 
       % reshape and extract red
       
       
       
  
        outpathRed  = fullfile(newpath,'redchannel.raw');
        outpathGreen= fullfile(newpath,'greenchannel.raw');
        
        if exist(outpathRed,'file')
            delete(outpathRed)
        end
        
        if exist(outpathGreen,'file')
            delete(outpathGreen)
        end 
        
        
        red_out = fopen(outpathRed,'w');
        fwrite(red_out, Red,'uint16');
        fclose(red_out);
        
        green_out = fopen(outpathGreen,'w');
        fwrite(green_out,Green,'uint16');     
        fclose(green_out);
        
    else
        warning('There were more channels than expected in %s',input.expname{i})
        
    end
   
    
end    
    fprintf('Elapsed time: %g. minutes \n', toc()/60);

end 
%%
    function filePreparation(input,i)
        % this is legacy code designed to create a new path for the extracted data
        % and move the appropriate
        % legacy code from Nik Francis 2016
        
       
        bb=strsplit(input.path,'\'); 
        for i = 1:length(input.expname)
             
            
            %Local path for data
             localpath = fullfile(input.path,input.expname{i}); 
            newpath = fullfile(input.savepath ,  bb{end}, ...
                input.expname{i}) ;
            mkdir(newpath)
            %ThorImage experiment params
            inpath = fullfile(input.path , input.expname{i}, ...
                'Experiment.xml');
            
            
            
            %Move relevant ThorImage files to local path
            Thor_path= fullfile(input.path ,...
                input.expname{i},...
                'Experiment.xml');
            
            copyfile(Thor_path,fullfile(newpath,'Experiment.xml'));
 
            % there should be either a thor timing file or H5 file depending on
            % using 2.1 or 3.0 respestively we move either one over here
            try
                Timing_path =fullfile(input.path,...
                    input.expname{i},...
                    'timing.txt');
                if ~exist( fullfile(newpath, 'timing.txt'),'file')
                    copyfile(Timing_path, fullfile(newpath,'timing.txt'));
                end
            catch
              H5 =   dir(fullfile(localpath, '**/*/Episode001.h5') );
                if ~isempty(H5)
                    h5file_path = fullfile(H5.folder,H5.name);
                    try
                    movefile(h5file_path,fullfile(newpath, 'Episode001.h5'));
                    catch
                    warning('H5 File Not moved')    
                    end 
                end
            end
        end
        
    end

