function AvgImage(input,channels)


if channels == 1 
for i = 1:size(input.datalist,1)

path = char(input.datalist(i,4)) ;   
opts = get_options_from_xml([path '/Experiment.xml']);
opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};
GreenChannel = memmapfile([path '/greenchannelregistered.raw'], ...
        'Format', opts.format, 'Repeat', opts.numframes);
greenChanImg=zeros(opts.dimX,opts.dimY);
GCref = GreenChannel.Data;

        for frame = 1:opts.numframes
             greenChanImg(:,:) = greenChanImg(:,:) + double(GCref(frame).channels);
   
%             greenChanImg(:,:,frame) = GreenChannel.Data(frame).channels;
%             redChanImg(:,:,frame) = GreenChannel.Data(frame).channels;
            %fprintf('loading frame %d/%d/n', frame, opts.numframes);
        end  
AvgGC = uint16(round(greenChanImg/opts.numframes));
path_green = [path '/AvgGC.raw'];
outfile{1} = fopen(path_green, 'w');
fwrite(outfile{1}, AvgGC, 'uint16');

end


elseif channels == 2

for i = 1:size(input.datalist,1)

path = char(input.datalist(i,4)) ;   
opts = get_options_from_xml([path '/Experiment.xml']);
opts.format = {'uint16', [opts.dimX, opts.dimY 1], 'channels'};


GreenChannel = memmapfile([path '/greenchannelregistered.raw'], ...
        'Format', opts.format, 'Repeat', opts.numframes);
RedChannel = memmapfile([path '/redchannelregistered.raw'], ...
        'Format', opts.format, 'Repeat', opts.numframes);

greenChanImg=zeros(opts.dimX,opts.dimY);
redChanImg=zeros(opts.dimX,opts.dimY);
GCref = GreenChannel.Data;
RCref = RedChannel.Data;
        for frame = 1:opts.numframes
             greenChanImg(:,:) = greenChanImg(:,:) + double(GCref(frame).channels);
             redChanImg(:,:) = redChanImg(:,:) + double(RCref(frame).channels);
%             greenChanImg(:,:,frame) = GreenChannel.Data(frame).channels;
%             redChanImg(:,:,frame) = GreenChannel.Data(frame).channels;
            %fprintf('loading frame %d/%d/n', frame, opts.numframes);
        end  

AvgGC = uint16(round(greenChanImg/opts.numframes));
AvgRC = uint16(round(redChanImg/opts.numframes));

path_green = [path '/AvgGC.raw'];
path_red = [path '/AvgRC.raw'];
outfile{1} = fopen(path_green, 'w');
outfile{2} = fopen(path_red, 'w');
fwrite(outfile{1}, AvgGC, 'uint16');
fwrite(outfile{2}, AvgRC, 'uint16');
end
end
else
end
