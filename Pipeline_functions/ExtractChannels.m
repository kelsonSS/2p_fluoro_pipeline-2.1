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

    disp(['Extracting ' char(input.expname) ])
    
    newpath = char(fullfile(input.savepath, input.expname));
    xmlpath = char(fullfile(newpath,'Experiment.xml'));
    opts = get_options_from_xml(xmlpath);
    
    if strcmp(opts.version(1:3),'4.0')
        img_name = 'Image_001_001.raw';
    else 
        img_name = 'Image_0001_0001.raw';
    end 
        
   
    out_path = char(fullfile(newpath,'greenchannel.raw'));
    
    % chcek to see if the file has already be extracted 
    if ~ isempty(dir(out_path))
        newfile = dir(out_path);
        oldfile = dir(out_path);         
        if newfile.bytes == oldfile.bytes 
            return
        end
    end 
    
    if opts.numchannels == 1
        
        old_path = fullfile(input.inpath,input.expname,img_name);
        try
        copyfile(old_path,out_path, 'f')
        catch
            return
        end
        
    elseif opts.numchannels == 2
        
        out_path = fullfile(newpath,'greenchannel');
        img = fopen(fullfile(input.inpath,input.expname,img_name));
       
        try
        Green = fread(img,'uint16=>uint16');
        catch 
            fclose(img);
            return
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
        
            
            
            
        else

            nchannels = 2 ;
            
          
            Red = Green(:,:,[2:2:end]);
            Green(:,:,[2:2:end]) = [];
        
        end 
       % reshape and extract red
       
       
       
  
        outpathRed  = fullfile(newpath,'redchannel.raw');
        outpathGreen= fullfile(newpath,'greenchannel.raw');
       
        
%         if exist(outpathRed,'file')
%             delete(outpathRed)
%         end
%         
%         if exist(outpathGreen,'file')
%             delete(outpathGreen)
%         end 
        
        
        red_out = fopen(outpathRed,'w');
        fwrite(red_out, Red,'uint16');
        fclose(red_out);
        
        green_out = fopen(outpathGreen,'w');
        fwrite(green_out,Green,'uint16');     
        fclose(green_out);
        
    else
        warning('There were more channels than expected in %s',input.expname)
        
    end
   
      
    fprintf('Elapsed time: %g. minutes \n', toc()/60);

end 
