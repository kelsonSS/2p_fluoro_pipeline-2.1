function ExtractGreenChannel(input)
%ExtractRedGreenChannels is a program that loads the raw image sequence
%from ThorImage and extracts just the red and green channels.
%Nikolas Fedit rancis 2016
tic();
%Load image sequences into memory
infiles=[];
FramesPerMovie=[];
bb=strsplit(input.path,'\'); %changed from strsep to strsplit (KA)
for i = 1:length(input.expname)
    
    %Local path for data
    newpath = [input.savepath '\'  bb{5}  '\' ...
        input.expname{i} '\'];
    mkdir(newpath)
    %ThorImage experiment params
    inpath = fullfile(input.path , input.expname{i}, ...
                      'Experiment.xml');
    opts = get_options_from_xml(inpath);
    opts.format = {'uint16', [opts.dimX, opts.dimY, 1], 'channels'};
    %Load image sequences
    infiles{i} = memmapfile([input.path '\' input.expname{i} ...
        '\Image_0001_0001.raw'], 'Format', opts.format, ...
        'Repeat', opts.numframes);
    FramesPerMovie = [FramesPerMovie; opts.numframes];
    %Move relevant ThorImage files to local path
    copyfile([input.path '\' input.expname{i} '\Experiment.xml'], ...
        [newpath 'Experiment.xml']);
    fprintf('copying %d of %d',i,length(input.expname) )
    try
    movefile([input.path '\' input.expname{i} '\timing.txt'], ...
        [newpath 'timing.txt']);
    catch
        if ~exist( [newpath 'Episode001.h5'],'file')
    movefile([input.path '\' input.expname{i} '\Episode001.h5'], ...
        [newpath 'Episode001.h5']);
        end 
    end 
end
%Extract red and green channels, and store them in the local path
for i = 1:length(input.expname)
    
    disp(['Extracting ' input.path '\' input.expname{i}])
    %Define raw sequence files in local path
    newpath = [input.savepath '\'  bb{5} '\' ...
        input.expname{i} '\'];
%     newpath_red = [newpath 'redchannel.raw'];
    newpath_green = [newpath 'greenchannel.raw'];
    
    % chcek to see if the file has already be extracted 
    if ~ isempty(dir(newpath_green))
        newfile = dir(newpath_green);
        oldfile = dir([input.path '\' input.expname{i}...
                      '\Image_0001_0001.raw']);         
        if newfile.bytes == oldfile.bytes 
            continue
        end
    end 
    
%   outfile{1} = fopen(newpath_red, 'w');
    outfile{1} = fopen(newpath_green, 'w');
    %Loop through frames and extract channels
    for frame = 1:FramesPerMovie(i)
        %redChanImg = infiles{i}.Data(frame).channels(:, :, opts.redChannel);
        greenChanImg = infiles{i}.Data(frame).channels(:, :, 1);
       % fwrite(outfile{1}, redChanImg, 'uint16');
        fwrite(outfile{1}, greenChanImg, 'uint16');
        %fprintf('frame %d/%d/n', frame, FramesPerMovie(i));
    end
    fclose(outfile{1});
    %fclose(outfile{2});
    %Copy Psignal data to local directory
%     copyfile([input.path '\' 'Psignal'], [input.savepath '\' bb{7} '\' bb{7} ...
%        '\Psignal'])
end
fprintf('Elapsed time: %g. minutes/n', toc()/60);