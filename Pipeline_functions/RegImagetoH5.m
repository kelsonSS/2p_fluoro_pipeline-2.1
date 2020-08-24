function RegImagetoH5(image_folder)
    % this function takes a registered raw image and packages it into and
    % H5file for use in suite2p. This functionality will eventually be
    % ported into the 2P-Pipeline this function can convert all old files
    % that were not converted 
    %
    % inputs
    % image_folder - The path to a folder that contains
    % 
    % Image - as greenchannelregistered.raw
    % xml - thorimage xml as Experiment.xml
    

 %% Load XML
  XMLPath  = fullfile(image_folder,'Experiment.xml');
  opts = get_options_from_xml(XMLPath);
 
 %% Load the input .raw file into a matrix
   FullPathIMG = fullfile(image_folder, 'greenchannel.raw');

fh = fopen(FullPathIMG); 
numFrames_total =  FindRawImgSize(fh,[opts.dimX opts.dimY]);
ItemsPerImage = opts.dimX* opts.dimY ;
IMG = fread(fh,inf,'uint16=>uint16','l');
IMG =  permute( reshape(IMG,opts.dimX, opts.dimY, []),[2 1 3]); 

%% create h5 file
H5Path = fullfile(image_folder,'greenChannel.h5');
if ~exist(H5Path,'file')
    h5create(H5Path,'/Images', [opts.dimY,opts.dimX,numFrames_total]);
end 
h5write(H5Path,'/Images',IMG)


