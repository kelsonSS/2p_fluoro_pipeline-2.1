function [smRoiBoundaries smNpBoundaries] = FindROIs(handles)
%Cell centers
xc = handles.selectedneurons.Data(:,1);
yc = handles.selectedneurons.Data(:,2);
%Contrast enhanced green channel
muIMG = handles.GreenContAdjfilteredadjMnIMG;
%Preallocate ROIs
pcimg = cell (length(xc) , 13 , 360 );
imgCrop = cell ( length(xc) , handles.cellcropdim, handles.cellcropdim );
imgCropNorm = cell (length(xc), handles.cellcropdim , handles.cellcropdim );
roiBoundaries = cell ( length(xc) , 360 , 3 );
smRoiBoundaries = cell ( length(xc) , 360 , 3 );
ROIOut = zeros ( 360 , 2 );
ROIIn =  zeros ( 360 , 2 );
%Preallocate Neuropil
npBoundaries = cell ( length(xc) , 1 , 3 );
smNpBoundaries = cell ( length(xc) , 1 , 3 );
NeuropilOut = zeros ( 360 , 2 );
NeuropilIn =  zeros ( 360 , 2 );
%Plot green channel
cla(handles.axes2)
imshow((handles.GreenContAdjfilteredadjMnIMG),[],'Parent',handles.axes2);
hold(handles.axes2,'on')
%Find boundary of neurons
for pp = 1:length(xc)
    %Crop the neuron
    xpt=xc(pp);
    ypt=yc(pp);
    imgCrop{pp} = imcrop(muIMG,[xpt-(handles.cellcropdim/2) ypt-(handles.cellcropdim/2) handles.cellcropdim handles.cellcropdim]);
    imgCropNorm{pp} = (imgCrop{pp}-min(imgCrop{pp}(:)))./(range(imgCrop{pp}(:)));
    %Convert the image to polar coordinates
    pcimg{pp} = imgpolarcoord(imgCropNorm{pp});
    %Find the ring inner edge. The ring edges are defined in the polor
    %coordinate space, where a ring will appear, ideally, as a thick line.
    %We look for the edge by finding where the top and bottom of the thick
    %line,
    RingPks = zeros(size(pcimg{pp},2),1);
    RingNeuropil = zeros(size(pcimg{pp},2),1);
    tmpNeuron = [];
    tmpNeuron = pcimg{pp};
    %Iterate over each of the 360 degrees
    for cc = 1:size(tmpNeuron,2)
        tmpThresh = [];
        thresh1 = [];
        pkTmp = [];
        %Find the peak for each degree
        pkTmp = find(diff(tmpNeuron(:,cc)) == max(diff(tmpNeuron(:,cc))));
        if ~isempty(pkTmp)
            %more than one pixel identified - grab the first one
            if length(pkTmp) > 1 
                if pkTmp(1) < handles.expectedNeuronRadiusPix
                    RingPks(cc) = pkTmp(1);
                else
                    RingPks(cc) = handles.expectedNeuronRadiusPix;
                end
            else
                if pkTmp < handles.expectedNeuronRadiusPix
                    RingPks(cc) = pkTmp;
                else
                    RingPks(cc) = handles.expectedNeuronRadiusPix;
                end
            end
            % if its'the first direction and no peaks are found
        elseif cc == 1 
            RingPks(cc) = handles.expectedNeuronRadiusPix; 
        else
            %If no peak was found within acceptable bounds, then use the last peak
            RingPks(cc) = RingPks(cc-1);
        end
        ROIOut(cc,:) = [ degtorad(cc)  RingPks(cc)+ceil(handles.ringthickness(1)) ];
        ROIIn(cc,:) = [ degtorad(cc)   RingPks(cc)];
        NeuropilIn(cc,:) = [degtorad(cc) ROIOut(cc,2)+handles.NPgap];
        NeuropilOut(cc,:) = [degtorad(cc) ROIOut(cc,2)+handles.NPthickness];
    end
    roiBoundaries{pp} = [ ROIIn(:,1) ROIIn(:,2)  ROIOut(:,2) ]; % [PolarCoords (0-2Pi)    InnerRing     OuterRing]
    smRoiBoundaries{pp} = [ ROIIn(:,1) smooth(ROIIn(:,2),handles.smoothfact)  smooth(ROIOut(:,2),handles.smoothfact) ]; % [PolarCoords (0-2pi)    InnerRing     OuterRing]
    npBoundaries{pp} = [ NeuropilIn(:,1) NeuropilIn(:,2) NeuropilOut(:,2) ];
    smNpBoundaries{pp} = [ NeuropilIn(:,1) smooth(NeuropilIn(:,2),handles.smoothfact) smooth(NeuropilOut(:,2),handles.smoothfact) ];
    plot(xpt,ypt,'r.','Parent',handles.axes2) %plot clicked centers
    plot(  xpt + smRoiBoundaries{pp}(:,3) .* (cos(smRoiBoundaries{pp}(:,1)))  ,  ypt + smRoiBoundaries{pp}(:,3) .* (sin(smRoiBoundaries{pp}(:,1))) , 'r','Parent',handles.axes2 ) %plot outer edges of rings
    plot(  xpt + smRoiBoundaries{pp}(:,2) .* (cos(smRoiBoundaries{pp}(:,1)))  ,  ypt + smRoiBoundaries{pp}(:,2) .* (sin(smRoiBoundaries{pp}(:,1))) , 'r','Parent',handles.axes2 ) %plot inner edges of rings
    plot(  xpt + smNpBoundaries{pp}(:,3) .* (cos(smNpBoundaries{pp}(:,1)))  ,  ypt + smNpBoundaries{pp}(:,3) .* (sin(smNpBoundaries{pp}(:,1))) , 'y' ,'Parent',handles.axes2) %plot outer edges of rings
    plot(  xpt + smNpBoundaries{pp}(:,2) .* (cos(smNpBoundaries{pp}(:,1)))  ,  ypt + smNpBoundaries{pp}(:,2) .* (sin(smNpBoundaries{pp}(:,1))) , 'y','Parent',handles.axes2 ) %plot inner edges of rings
end