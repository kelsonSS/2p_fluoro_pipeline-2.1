function varargout = CellDefinitionGUI(varargin)
% if ~ismac
%     javax.swing.UIManager.setLookAndFeel('com.sun.java.swing.plaf.windows.WindowsLookAndFeel')
% else
%     javax.swing.UIManager.setLookAndFeel('com.sun.java.swing.plaf.com.apple.laf.AquaLookAndFeel')
% end
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CellDefinitionGUI_OpeningFcn, ...
    'gui_OutputFcn',  @CellDefinitionGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function CellDefinitionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
set(hObject,'Position',[50,24,4000,50])
handles.output = hObject;
handles.expectedNeuronDiamMicrons = varargin{1}.expectedNeuronDiamMicrons;
%Load parameters from input variable
handles.cellcropdim = varargin{1}.cellcropdim;
handles.ringthickness = varargin{1}.ringthickness;
handles.NPthickness = varargin{1}.NPthickness;
handles.NPgap = varargin{1}.NPgap;
handles.smoothfact = varargin{1}.smoothfact;
cla(handles.axes1)
cla(handles.axes2)
guidata(hObject,handles);
function varargout = CellDefinitionGUI_OutputFcn(hObject, eventdata, handles)
guidata(hObject,handles);
function loadimage_Callback(hObject, eventdata, handles)
%Select raw image to segement
if isfield(handles,'filteredadjMnIMG')
    handles=rmfield(handles,'filteredadjMnIMG');
end
cla(handles.axes1)
title('')
text(.1,.5,'LOADING IMAGES...','fontsize',20)
[fname, pthname] = uigetfile('\\vault3\data\kelson\analyzed\*.raw','Select a RED channel DataFile');

    disp('User selected Cancel');
    handles.pthname = 0;
    guidata(hObject,handles);
    fullpathIMG = [pthname fname];
    RfullpathIMG = [pthname fname];
    cla(handles.axes1)
    

    [Gfname, Gpthname] = uigetfile([pthname '\*.raw'], ...
    'Select a GREEN channel DataFile');
    fullpathIMG = [pthname fname];
    GfullpathIMG = [Gpthname Gfname];
    
    handles.pthname = {pthname};
    bb=strsplit(pthname,'\');
    bb = bb{end-1};
    handles.filename = [bb];
    guidata(hObject,handles);

if iscellstr({pthname})
    %Read associated XML file for image dimensions
    fullpathXML = [pthname 'Experiment.xml'];
    opts = get_options_from_xml(fullpathXML);
    numChannels = 1;
    numTimePts = opts.numframes;
   % numImages = min([1 numChannels * numTimePts]);
    numImages = min([2000 numChannels * numTimePts]);
    dimX = opts.dimX;
    dimY = opts.dimY;
    %Read RAW file into workspace
    figure(handles.figure1)
    %Process Red image
    fh = fopen(RfullpathIMG);

    array = fread(fh,[dimX*dimY*numImages],'uint16=>uint16');
    
    numImagesActual = length(array) / dimX /dimY;
    
    numImages = min(numImages,numImagesActual);
    IMG = reshape(array, [dimX, dimY, numImages]);
    red_rotate_flag = questdlg('Rotate Red IMG?') 
    fclose(fh);
    IMG = permute(IMG, [2 1 3]);

  
    
    mnIMG = squeeze(mean(IMG(:,:,1:end),3));
    adjMnIMG = (mnIMG - min(mnIMG(:))) ./ (range(mnIMG(:)));
     
      
      
    %Process green image
    fh = fopen(GfullpathIMG);
    numImages = min([2000 numChannels * numTimePts]);
    
    array = fread(fh,[dimX*dimY*numImages],'uint16=>uint16');
    numImagesActual = length(array) / dimX /dimY;
    numImages = min(numImages,numImagesActual); 
    
    
    IMG = reshape(array,[dimX, dimY,numImages]);
    IMG = permute(IMG, [2 1 3]);
    fclose(fh);
    mnIMG = squeeze(mean(double(IMG(:,:,:)),3));
    GadjMnIMG = (mnIMG - min(mnIMG(:))) ./ (range(mnIMG(:)));
    [avgRedfilteredImage filteredImages] = hmfilter(adjMnIMG);
    [avgGreenfilteredImage filteredImages] = hmfilter(GadjMnIMG);
    handles.RedfilteredadjMnIMG=avgRedfilteredImage;
    handles.GreenfilteredadjMnIMG=avgGreenfilteredImage;
    try
    fNameString = ['AvgGreenImage.mat'];
    idx = strfind(handles.pthname{1},'\');
    pthname = handles.pthname{1}(1:idx(end-1));
    DestPath=fullfile(pthname,handles.filename,fNameString);
    save(DestPath,'avgGreenfilteredImage');
    catch
        warning('unable to save green image' ) 
    end 
    
    
    %Plot images
    figure(handles.figure1)
    handles.RedContAdjfilteredadjMnIMG = handles.RedfilteredadjMnIMG;
    handles.GreenContAdjfilteredadjMnIMG = handles.GreenfilteredadjMnIMG;
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    set(handles.contrast_red,'Max',0.49);
    set(handles.contrast_red,'Value',0.49);
    title(handles.axes1, [handles.filename ': Red Channel'])
    cla(handles.axes2)
    imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
    set(handles.contrast_green,'Max',0.49);
    set(handles.contrast_green,'Value',0.49);
    title(handles.axes2, [handles.filename ': Green Channel'])
    set(handles.figure1,'toolbar','figure');
    set(handles.figure1,'menubar','figure');
    
%    if exist(fullfile(pthname,handles.filename,'AvgImageSmooth.tif'),'file')
%        figure
%        fh = Tiff(fullfile(pthname,handles.filename,'AvgImageSmooth.tif'),'r')
%        img = fh.read();
%        img = permute(img,[2 ,1,3]);
%        fh.close()
%        load(fullfile(pthname,handles.filename,'smooth_img.mat'))
%        imshow(permute(smooth_img,[2,1,3]) ,[],'Parent',handles.axes2);
%        smooth_img = (smooth_img - min(smooth_img(:))) ./ (range(smooth_img(:)));
%     end 
%    handles.GreenContAdjfilteredadjMnIMG = smooth_img; 
%    handles.GreenfilteredadjMnIMG = smooth_img;
else 
    cla(handles.axes1)
end
guidata(hObject,handles)
function contrast_red_Callback(hObject, eventdata, handles)
handles.numonoff.Value=1;
contval = get(handles.contrast_red,'Value');
if isfield(handles,'RedfilteredadjMnIMG')
    handles.RedContAdjfilteredadjMnIMG = imadjust(handles.RedfilteredadjMnIMG, ...
        [0; contval+eps],[]);
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    title(handles.axes1, [handles.filename ': Red Channel'])
    hold(handles.axes1,'on')
    try
        xy= handles.selectedneurons.Data;
        xc = xy(:,1);
        yc = xy(:,2);
        ptsIdx = [[1:length(xc)]' xc yc];
        title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes1)
        for nn = 1:length(xc)
            hold on
            neuronNumb = num2str(ptsIdx(nn));
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes1)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.RedContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
        end
        guidata(hObject,handles);
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ...
            ' Neuron(s) Selected'],'Parent',handles.axes1)
    end
end
guidata(hObject, handles);
function contrast_green_Callback(hObject, eventdata, handles)
handles.numonoff.Value = 1;
contval = get(handles.contrast_green,'Value');
if isfield(handles,'GreenfilteredadjMnIMG')
    handles.GreenContAdjfilteredadjMnIMG = imadjust(handles.GreenfilteredadjMnIMG, ...
        [0; contval+eps],[]);
    cla(handles.axes2)
    imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
    title(handles.axes2, [handles.filename ': Green Channel'])
    hold(handles.axes2,'on')
    try
        xy= handles.selectedneurons.Data;
        xc = xy(:,1);
        yc = xy(:,2);
        ptsIdx = [[1:length(xc)]' xc yc];
        title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes2)
        for nn = 1:length(xc)
            hold on
            neuronNumb = num2str(ptsIdx(nn));
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes2)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes2)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.GreenContAdjfilteredadjMnIMG,1)/32,'Parent', ...
                handles.axes2)
        end
        guidata(hObject,handles);
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ...
            ' Neuron(s) Selected'],'Parent',handles.axes2)
    end
end
guidata(hObject, handles);
function selectneurons_Callback(hObject, eventdata, handles)
%handles.numonoff.Value=1;
%Select neuron centers
title(['0 Neurons Selected'],'Parent',handles.axes1)
title(['0 Neurons Selected'],'Parent',handles.axes2)
figure(handles.figure1)
if isfield(handles,'RedfilteredadjMnIMG') && iscellstr(handles.pthname)
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    hold(handles.axes1,'on')
    cla(handles.axes2)
    imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
    hold(handles.axes2,'on')
    try
        if handles.RG.Value
            [xc, yc] = getpts(handles.axes2);
        else
            [xc, yc] = getpts(handles.axes1);
        end
        Oframe = sum(xc<0)+ ...
            sum(xc>size(handles.RedContAdjfilteredadjMnIMG,1))+ ...
            sum(yc<0)+sum(yc>size(handles.RedContAdjfilteredadjMnIMG,2));
        for i = 1:length(xc)
            if sum(xc(i)<0)+ ...
                sum(xc(i)>size(handles.RedContAdjfilteredadjMnIMG,1))+ ...
                sum(yc(i)<0)+sum(yc(i)>size(handles.RedContAdjfilteredadjMnIMG,2));
                xc(i)=[];
                yc(i)=[];
            end
        end
        if Oframe > 0
            warndlg('Selection(s) out of frame! Outliers removed...')
            if isempty(xc)
                return
            end
        end
        ptsIdx = [[1:length(xc)]' xc yc];
        title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes1)
        title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes2)
        for nn = 1:length(xc)
            hold on
            neuronNumb = num2str(ptsIdx(nn));
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes2)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes2)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.RedfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.GreenfilteredadjMnIMG,1)/32,'Parent',handles.axes2)
        end
        handles.selectedneurons.Data = [xc yc];
        guidata(hObject,handles);
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ...
            ' Neuron(s) Selected'],'Parent',handles.axes1)
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ...
            ' Neuron(s) Selected'],'Parent',handles.axes2)
    end
end
guidata(hObject, handles);
function deleteselections_Callback(hObject, eventdata, handles)
try
    handles.selectedneurons.Data(handles.deleteneuron(:,1),:)=[];
end
figure(handles.figure1)
if isfield(handles,'RedContAdjfilteredadjMnIMG') && iscellstr(handles.pthname)
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    hold(handles.axes1,'on')
    cla(handles.axes2)
    imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
    hold(handles.axes2,'on')
    rois_Callback(hObject,eventdata, handles)
    try
        xy= handles.selectedneurons.Data;
        xc = xy(:,1);
        yc = xy(:,2);
        ptsIdx = [[1:length(xc)]' xc yc];
        title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes1)
        title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes2)
        for nn = 1:length(xc)
            hold on
            neuronNumb = num2str(ptsIdx(nn));
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
            'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes2)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes2)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.RedContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
            plot(xc(nn), yc(nn),'r.','markersize', ...
            size(handles.GreenContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes2)
        end
        guidata(hObject,handles);
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ...
            ' Neuron(s) Selected'],'Parent',handles.axes1)
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ...
            ' Neuron(s) Selected'],'Parent',handles.axes2)
    end
end
guidata(hObject,handles);
function addneurons_Callback(hObject, eventdata, handles)
try
    if handles.RG.Value
        [xc, yc] = getpts(handles.axes2);
    else
        [xc, yc] = getpts(handles.axes1);
    end
    Oframe = sum(xc<0)+ ...
        sum(xc>size(handles.RedContAdjfilteredadjMnIMG,1))+ ...
        sum(yc<0)+sum(yc>size(handles.RedContAdjfilteredadjMnIMG,2));
    for i = 1:length(xc)
        if sum(xc(i)<0)+ ...
                sum(xc(i)>size(handles.RedContAdjfilteredadjMnIMG,1))+ ...
                sum(yc(i)<0)+sum(yc(i)>size(handles.RedContAdjfilteredadjMnIMG,2));
            xc(i)=[];
            yc(i)=[];
        end
    end
    if Oframe > 0
        warndlg('Selection(s) out of frame! Outliers removed...')
        if isempty(xc)
            return
        end
    end
    xc = [handles.selectedneurons.Data(:,1); xc];
    yc = [handles.selectedneurons.Data(:,2); yc];
    ptsIdx = [[1:length(xc)]' [xc yc]];
    handles.selectedneurons.Data = [xc yc];
    
    title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes1)
    title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes2)
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    hold(handles.axes1,'on')
    cla(handles.axes2)
    imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
    hold(handles.axes2,'on')
%     rois_Callback(hObject,eventdata, handles)
    guidata(hObject,handles)
    numonoff_Callback(hObject, eventdata, handles)
   

    title([handles.filename ': ' num2str(size(ptsIdx,1)) ' Neuron(s) Selected'], ...
        'Parent',handles.axes1)
    title([handles.filename ': ' num2str(size(ptsIdx,1)) ' Neuron(s) Selected'], ...
        'Parent',handles.axes2)
end
function selectedneurons_CellSelectionCallback(hObject, eventdata, handles)
handles.deleteneuron = eventdata.Indices;
guidata(hObject,handles);
function numonoff_Callback(hObject, eventdata, handles)
    
    xc = handles.selectedneurons.Data(:,1);
    yc = handles.selectedneurons.Data(:,2);
    ptsIdx = [[1:length(xc)]' [xc yc]];
    title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes1)
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    hold(handles.axes1,'on')
    title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes2)
    cla(handles.axes2)
    imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
    hold(handles.axes2,'on')
    rois_Callback(hObject,eventdata, handles)
    
    
    
    
    
    if ~handles.numonoff.Value
    for nn = 1:length(xc)
        hold on
        plot(xc(nn), yc(nn),'r.','markersize', ...
            size(handles.RedContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
        plot(xc(nn), yc(nn),'r.','markersize', ...
            size(handles.GreenContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes2)
    end
    guidata(hObject,handles);
    title([handles.filename ': ' num2str(size(ptsIdx,1)) ' Neuron(s) Selected'], ...
        'Parent',handles.axes1)
    title([handles.filename ': ' num2str(size(ptsIdx,1)) ' Neuron(s) Selected'], ...
        'Parent',handles.axes2)
    else
        
        for nn = 1:length(xc)
            hold on
            neuronNumb = num2str(ptsIdx(nn));
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes1)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',16,'Color','b', ...
                'fontweight','bold','Parent',handles.axes2)
            text(xc(nn)+10, yc(nn)+10,neuronNumb,'FontSize',14,'Color','c', ...
                'fontweight','bold','Parent',handles.axes2)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.RedContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
            plot(xc(nn), yc(nn),'r.','markersize', ...
                size(handles.GreenContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes2)
        end
        guidata(hObject,handles);
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ' Neuron(s) Selected'], ...
            'Parent',handles.axes1)
        title([handles.filename ': ' num2str(size(ptsIdx,1)) ' Neuron(s) Selected'], ...
            'Parent',handles.axes2)
end
function loadselection_Callback(hObject, eventdata, handles)
fNameString = ['CellDefinitions.mat'];
idx = strfind(handles.pthname{1},'\');
pthname = handles.pthname{1}(1:idx(end-1));
DestPath=fullfile(pthname, handles.filename,fNameString);
load(DestPath);
handles.selectedneurons.Data = ptsIdx(:,2:3);
function saveselection_Callback(hObject, eventdata, handles)
try
    handles.expectedNeuronRadiusPix=round(handles.expectedNeuronDiamMicrons/handles.micronsPerPixelX)/2;
    [smRoiBoundaries smNpBoundaries] = FindROIs(handles);
    handles.smRoiBoundaries = smRoiBoundaries;  
    handles.smNpBoundaries = smNpBoundaries;
    
    ptsIdx = [[1:size(handles.selectedneurons.Data,1)]' handles.selectedneurons.Data];
    fNameString = ['CellDefinitions.mat'];
    idx = strfind(handles.pthname{1},'\');
    pthname = handles.pthname{1}(1:idx(end-1));
    DestPath=fullfile(pthname,handles.filename,fNameString);
    Cdef=[];
    Cdef.ptsIdx = ptsIdx;
    Cdef.contrast = get(handles.contrast_green,'Value');
    Cdef.expectedNeuronDiamMicrons = handles.expectedNeuronDiamMicrons;
    Cdef.expectedNeuronRadiusPix = handles.expectedNeuronRadiusPix;
    Cdef.cellcropdim = handles.cellcropdim;
    Cdef.ringthickness = handles.ringthickness;
    Cdef.NPgap = handles.NPgap;
    Cdef.smoothfact = handles.smoothfact;
    Cdef.smRoiBoundaries = handles.smRoiBoundaries;
    Cdef.smNpBoundaries = handles.smNpBoundaries;
    Cdef.dimXmicrons = handles.dimXmicrons;
    Cdef.micronsPerPixelX = handles.micronsPerPixelX;
    Cdef.dimX = handles.dimX;
    Cdef.ringthickness = handles.ringthickness;
    Cdef.NPgap = handles.NPgap;
    Cdef.smoothfact = handles.smoothfact;
    save(DestPath, '-struct', 'Cdef');
catch
    warndlg('NO NEURONS SELECTED')

end
function [avgfilteredImage filteredImages] = hmfilter(I)
filteredImages=zeros(size(I));

for i = 1:size(I,3)
    I_temp = I(:,:,i);
    I_temp = im2double(I_temp);
    I_temp = log(1 + I_temp);
    M = 2*size(I_temp,1);
    N = 2*size(I_temp,2);
    sigma = 1.5; 
    [X, Y] = meshgrid(1:N,1:M);
    centerX = ceil(N/2);
    centerY = ceil(M/2);
    gaussianNumerator = (X - centerX).^2 + (Y - centerY).^2;
    H = exp(-gaussianNumerator./(2*sigma.^2));
    H=1-H;
    H = fftshift(H);
    Ir = padarray(I_temp,[ceil(size(I_temp,1)/2) ceil(size(I_temp,2)/2)],'symmetric');
    If = fft2(Ir, M, N);
    Iout = real(ifft2(H.*If));
    Iout = Iout(ceil(size(I_temp,1)/2)+1:size(Iout,1)-ceil(size(I_temp,1)/2), ...
        ceil(size(I_temp,2)/2)+1:(size(Iout,2)-ceil(size(I_temp,2)/2)));
    if i == 1 
     filteredImages(:,:) = exp(Iout) - 1;
    else 
    filteredImages(:,:,i) = exp(Iout) - 1;
    end 
end
avgfilteredImage = squeeze(mean(filteredImages,3));
function rois_Callback(hObject, eventdata, handles)
%Read associated XML file for image dimensions
fullpathXML = [handles.pthname{1} 'Experiment.xml'];
opts = get_options_from_xml (fullpathXML);
handles.dimX = opts.dimX;
handles.dimXmicrons = opts.dimXmicrons;
handles.micronsPerPixelX = opts.micronsPerPixel;

handles.expectedNeuronRadiusPix=round(handles.expectedNeuronDiamMicrons/handles.micronsPerPixelX)/2;
[smRoiBoundaries smNpBoundaries] = FindROIs(handles);
handles.smRoiBoundaries = smRoiBoundaries;
handles.smNpBoundaries = smNpBoundaries;
guidata(hObject,handles);
function info_Callback(hObject, eventdata, handles)
helpdlg([{['Cell Definition v1.0. Written by Dan Winkowski and Nikolas Francis.']};{' '};
    {['Cell Definition is ',...
    'used to find the coordinates of neurons in data collected from 2 photon imaging. ' ... 
    'For single-channel recordings, the user loads the same channel twice, or for red-labeled ' ...
    'mice, the user loads the red and green channels. The user then selects/adds/',...
    'deletes cell locations while dynamically adjusting contrast. The cell definition file is ' ...
    'saved in the same directory as the loaded images.']}],'Cell Definition HELP');
function contrast_red_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function contrast_green_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function RG_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(hObject, 'BackgroundColor',[0 1 0])
else
    set(hObject, 'BackgroundColor',[1 0 0])
end


function MovePoints_Callback(hObject, eventdata, handles)
 % plot points in preperation for KeyPressFcn Callback
xc = handles.selectedneurons.Data(:,1);
yc = handles.selectedneurons.Data(:,2);
    
    cla(handles.axes1)
    imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
    hold(handles.axes1,'on')

cla(handles.axes2)
imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
hold(handles.axes2,'on')
plot(xc, yc,'r.','markersize', ...
    size(handles.RedContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
plot(xc, yc,'r.','markersize', ...
    size(handles.GreenContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes2)
 
    
    


% --- Executes on key press with focus on pushbutton17 and none of its controls.
function MovePoints_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

xc = handles.selectedneurons.Data(:,1);
yc = handles.selectedneurons.Data(:,2);
ptsIdx = [[1:length(xc)]' [xc yc]];
switch eventdata.Key
        case 'leftarrow' 
            xc = xc-1;
        case 'rightarrow'
            xc = xc+1;
        case 'uparrow'
            yc  = yc-1;
        case 'downarrow'
           yc = yc+1;

end 
            
   % plotting
title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes1)
cla(handles.axes1)
imshow((handles.RedContAdjfilteredadjMnIMG),[],'Parent',handles.axes1);
hold(handles.axes1,'on')
title([num2str(size(ptsIdx,1)) ' Neurons Selected'],'Parent',handles.axes2)
cla(handles.axes2)
imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
hold(handles.axes2,'on')
plot(xc, yc,'r.','markersize', ...
    size(handles.RedContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes1)
plot(xc, yc,'r.','markersize', ...
    size(handles.GreenContAdjfilteredadjMnIMG,1)/32,'Parent',handles.axes2)
 

% cleanup
handles.selectedneurons.Data(:,1) = xc;
handles.selectedneurons.Data(:,2) = yc;

guidata(hObject,handles)


% --- Executes on button press in SaveNeuronCount.
function SaveNeuronCount_Callback(hObject, eventdata, handles)
% hObject    handle to SaveNeuronCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
  %% get number of neurons currently selected 
n_neurons = length(handles.selectedneurons.Data(:,1));

%% saving
fString = 'RedNeuronNumber.mat'
idx = strfind(handles.pthname{1},'\');
pthname = handles.pthname{1}(1:idx(end-1));
DestPath = fullfile(pthname,handles.filename, fString)

save(DestPath,'n_neurons')


% --- Executes on button press in Overlay_Smooth.
function Overlay_Smooth_Callback(hObject, eventdata, handles)
% hObject    handle to Overlay_Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   smooth_file = fullfile(handles.pthname{1},'AvgImageSmooth.tif');
      if exist(smooth_file,'file')
       %figure
        fh = Tiff(smooth_file,'r')
        img = fh.read();
        img = permute(img,[2 ,1,3]);
        fh.close()
        load(fullfile(handles.pthname{1},'smooth_img.mat'))
         smooth_img = (smooth_img - min(smooth_img(:))) ./ (range(smooth_img(:)));
        imshow(permute(smooth_img,[2,1,3]) ,[],'Parent',handles.axes2);
      end  

guidata(hObject,handles)
     

