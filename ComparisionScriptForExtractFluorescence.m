% This script proves the equality of the outputs of the previous image
% extraction method and my vectorized method
% previous method
img_path = '\\VAULT2\Vault2Data\Kelson\b248\2018-08-27\fra\Image_0001_0001.raw';


% previous method
m = memmapfile(img_path,'Format', opts.format,'Repeat', 12150); 
tic();for ii = 1:12150
img1(:,:,ii) = m.Data(ii).channels(:,:,1);
end
time1 = toc();
% my method
img = fopen('\\VAULT2\Vault2Data\Kelson\b248\2018-08-27\fra\Image_0001_0001.raw','r');

tic();img2 = fread(img,'uint16=>uint16');img2=reshape(img1,512,512,[]);
time2 = toc();
fclose(img);


isequal(img1,img2)

 