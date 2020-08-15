function exptVars = xmlVersionCheck(XML,numImagesOverwrite)

% Check ThorImage version for correct pointers -- could also use a strfind
% but some of the str names changed with newer versions therefore you have
% to check every time a new version comes out anyway...
if strcmp(XML.Children(10).Attributes.Value,'4.0.2019.8191')
    fprintf('Experiment used ThorImage Version %s\n',XML.Children(10).Attributes.Value)
    exptVars.dimX = str2double(XML.Children(26).Attributes(37).Value); % in pixels
    exptVars.dimY = str2double(XML.Children(26).Attributes(38).Value);
    exptVars.numImages = str2double(XML.Children(42).Attributes(7).Value);
    exptVars.micronsPerPixel = str2double(XML.Children(26).Attributes(36).Value);
    exptVars.frameRate = round(str2double(XML.Children(26).Attributes(24).Value)); % Not 30?
    exptVars.flybackFrames = str2double(XML.Children(42).Attributes(4).Value);
    exptVars.stepSizeUM = str2double(XML.Children(18).Attributes(5).Value);
    exptVars.numZSteps = str2double(XML.Children(18).Attributes(6).Value);
    exptVars.totalZplanes = exptVars.numZSteps+exptVars.flybackFrames;
    if exptVars.numZSteps==1
        exptVars.totalZplanes = 1;
    end
elseif strcmp(XML.Children(10).Attributes.Value,'3.1.2017.10021')
    fprintf('Experiment used ThorImage Version %s\n',XML.Children(10).Attributes.Value)
    exptVars.dimX = str2double(XML.Children(26).Attributes(31).Value); % in pixels
    exptVars.dimY = str2double(XML.Children(26).Attributes(32).Value);
    exptVars.numImages = str2double(XML.Children(42).Attributes(5).Value);
    exptVars.micronsPerPixel = str2double(XML.Children(26).Attributes(30).Value);
    exptVars.frameRate = round(str2double(XML.Children(26).Attributes(19).Value));
    zScanEnabled = str2double(XML.Children(42).Attributes(13).Value);
    if zScanEnabled
        exptVars.flybackFrames = str2double(XML.Children(42).Attributes(4).Value);
        exptVars.stepSizeUM = str2double(XML.Children(18).Attributes(5).Value);
        exptVars.numZSteps = str2double(XML.Children(18).Attributes(6).Value);
        exptVars.totalZplanes = exptVars.numZSteps+exptVars.flybackFrames;
    else
        exptVars.flybackFrames = 0;
        exptVars.stepSizeUM = 0;
        exptVars.numZSteps = 1;
        exptVars.totalZplanes = 1;
    end
elseif strcmp(XML.Children(10).Attributes.Value,'3.0.2016.10131')
    fprintf('Experiment used ThorImage Version %s\n',XML.Children(10).Attributes.Value)
    exptVars.dimX = str2double(XML.Children(26).Attributes(23).Value); % in pixels
    exptVars.dimY = str2double(XML.Children(26).Attributes(24).Value);
    exptVars.numImages = str2double(XML.Children(34).Attributes(5).Value); % FIX ***************
    exptVars.micronsPerPixel = str2double(XML.Children(26).Attributes(22).Value);
    exptVars.frameRate = round(str2double(XML.Children(26).Attributes(11).Value));
    exptVars.flybackFrames = str2double(XML.Children(34).Attributes(4).Value);
    exptVars.stepSizeUM = str2double(XML.Children(18).Attributes(3).Value);
    exptVars.numZSteps = str2double(XML.Children(18).Attributes(4).Value);
    exptVars.totalZplanes = exptVars.numZSteps+exptVars.flybackFrames;
    if exptVars.numZSteps==1
        exptVars.totalZplanes = 1;
    end
elseif strcmp(XML.Children(10).Attributes.Value,'2.1.2014.9102') % Not using old ThorImage moving forward
    fprintf('Experiment used ThorImage Version %s\n',XML.Children(10).Attributes.Value)
    exptVars.dimX = str2double(XML.Children(28).Attributes(16).Value); % in pixels
    exptVars.dimY = str2double(XML.Children(28).Attributes(16).Value);
    exptVars.dimXmicrons = str2double(XML.Children(22).Children(2).Children(2).Attributes(2).Value) *1000;
    exptVars.micronsPerPixel = exptVars.dimXmicrons/exptVars.dimX;
    exptVars.frameRate = 30; % FIND THIS XML POINTER
    exptVars.numImages = str2double(XML.Children(36).Attributes(5).Value);
    exptVars.flybackFrames = 0;
    exptVars.numZSteps = 1;
    exptVars.totalZplanes = exptVars.numZSteps+exptVars.flybackFrames;
else
    fprintf('Experiment used unfamiliar ThorImage Version %s, please add to code\n',XML.Children(10).Attributes.Value)
end

if exist('numImagesOverwrite','var') && ~isempty(numImagesOverwrite)
    exptVars.numImages = numImagesOverwrite;
end

%% Version-invariant parameters
exptVars.numVolumes = round(exptVars.numImages/exptVars.totalZplanes);
numBytes = 2; % assumes 16-bit
exptVars.segmentSize = (exptVars.dimX * exptVars.dimY * numBytes); %size of image in bytes
        
end