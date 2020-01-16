function ExtractRedGreenChannels(input)
%ExtractRedGreenChannels is a program that loads the raw image sequence
%from ThorImage and extracts just the red and green channels.
%Nikolas Francis 2016
tic();
%Load image sequences into memory
infiles=[];
FramesPerMovie=[];
bb=strsplit(input.path,'/'); %changed from strsep to strsplit (KA)
for i = 1:length(input.expname)
    %Local path for data
    newpath = [input.savepath '/' bb{8} '/' bb{9} '/' ...
        input.expname{i} '/'];
    mkdir(newpath)
    %ThorImage experiment params
    opts = get_options_from_xml([input.path '/' input.expname{i} ...
        '/Experiment.xml']);
    %Load image sequences
    infiles{i} = memmapfile([input.path '/' input.expname{i} ...
        '/Image_0001_0001.raw'], 'Format', opts.format, ...
        'Repeat', opts.numframes);
    FramesPerMovie = [FramesPerMovie; opts.numframes];
    %Copy relevant ThorImage files to local path
%     copyfile([input.path '/' input.expname{i} '/Experiment.xml'], ...
%         [newpath 'Experiment.xml']);
%     copyfile([input.path '/' input.expname{i} '/timing.txt'], ...
%         [newpath 'timing.txt']);
end
%Extract red and green channels, and store them in the local path
for i = 1:length(input.expname)
    disp(['Extracting ' input.path '/' input.expname{i}])
    %Define raw sequence files in local path
    newpath = [input.savepath '/' bb{8} '/' bb{9} '/' ...
        input.expname{i} '/'];
    newpath_red = [newpath 'redchannel.raw'];
    newpath_green = [newpath 'greenchannel.raw'];
    outfile{1} = fopen(newpath_red, 'w');
    outfile{2} = fopen(newpath_green, 'w');
    %Loop through frames and extract channels
    for frame = 1:FramesPerMovie(i),
        redChanImg = infiles{i}.Data(frame).channels(:, :, opts.redChannel);
        greenChanImg = infiles{i}.Data(frame).channels(:, :, opts.greenChannel);
        fwrite(outfile{1}, redChanImg, 'uint16');
        fwrite(outfile{2}, greenChanImg, 'uint16');
        %fprintf('frame %d/%d/n', frame, FramesPerMovie(i));
    end
    fclose(outfile{1});
    fclose(outfile{2});
    %Copy Psignal data to local directory
%     copyfile([input.path '/' 'Psignal'], [input.savepath '/' bb{8} '/' bb{8} ...
%        '/Psignal'])
end
fprintf('Elapsed time: %g. minutes/n', toc()/60);