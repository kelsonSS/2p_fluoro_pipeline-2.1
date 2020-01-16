function I = Tifload(deltaFfile)


I=[];
framecount=0;
InfoImage=imfinfo(deltaFfile);
NumberImages=length(InfoImage);
TifLink = Tiff(deltaFfile, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    framecount=framecount+1;
    I(:,:,framecount)=imresize(rot90(TifLink.read()),1);
    disp(['Loaded frame ' num2str(framecount)])
end
TifLink.close();