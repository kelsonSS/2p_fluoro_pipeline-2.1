function opts = get_options_from_xml(xmlpath,input)
%Parse the 'Experiment.xml' file at xmlpath and construct options struct opts from
%its contents.

opts = get_default_options;
xml =  xml2struct(xmlpath);
xml = xml.ThorImageExperiment;

opts.numframes = str2double(xml.Timelapse.Attributes.timepoints);
opts.stopFrame = opts.numframes;
opts.sections = 1;               
opts.dimX =  str2double(xml.LSM.Attributes.pixelX);
opts.dimY =  str2double(xml.LSM.Attributes.pixelY);
opts.TriggerMode=xml.Timelapse.Attributes.triggerMode; % 0 for finite, 2 for stimulus
opts.Streaming = xml.Streaming.Attributes.triggerMode; % 4 for fininte, 1 for stimulus
opts.StimTrigger = xml.Streaming.Attributes.stimulusTriggering;% 0 for finite 1 for stimulus?   
try
opts.micronsPerPixel = str2double(xml.LSM.Attributes.pixelSizeUM);
opts.dimXmicrons = str2double(xml.LSM.Attributes.widthUM);
catch
opts.dimXmicrons = str2double(xml.Sample.Wells.SubImages.Attributes.subOffsetXMM)*1000;
opts.micronsPerPixel = opts.dimXmicrons/opts.dimX;
end
opts.numchannels = length(xml.Wavelengths.Wavelength);
opts.format ={'uint16',...
             [opts.dimX,opts.dimY,opts.numchannels],...
             'channels'};
         
 opts.version = xml.Software.Attributes.version;


%The green channel will be either the first or second
%channel. Whether or not the channel is enabled depends on
%whether the corresponding gain is positive. This section of
%codes determines which channel is the green channel.
pmt = xml.PMT.Attributes;
gainA = pmt.gainA;
gainB = pmt.gainB;
assert(xor(strcmp(gainA, '0'), strcmp(gainB, '0')));
if ~strcmp(gainA, '0')
    opts.greenChannel = 1;
elseif ~strcmp(gainB, '0')
    opts.greenChannel = 2;
end
%As above, but for the red channel.
gainC = pmt.gainC;
gainD = pmt.gainD;
%assert(xor(strcmp(gainC, '0'), strcmp(gainD, '0')));
if ~strcmp(gainC, '0')
    opts.redChannel = 2;
elseif ~strcmp(gainD, '0')
    opts.redChannel = 4;
end



function opts = get_default_options
%Default input options for loading raw image sequences from ThorImage.
opts = struct;
%Image dimensions
opts.dimX = 512;
opts.dimY = 512;
%Channel index
opts.redChannel = 2;
opts.greenChannel = 1;
opts.templateFrame = 1;
opts.numframes = 5000; %Needs to be changed based on length of seq
%File formatting options
opts.format = {'uint16', [opts.dimX, opts.dimY, 2], 'channels'};
opts.startFrame = 1;
opts.stopFrame = opts.numframes;
opts.echo = true;
opts.profile = true;