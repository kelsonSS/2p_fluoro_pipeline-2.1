fclose all;
img = fopen('\\vault2\Vault2Data\Kelson\B264\18-11-23\OFC_silence\Image_0001_0001.raw');
Green = fread(img,'uint16=>uint16');
Green = reshape(Green, 512,512,[]);
%Red = Green(:,:,2:2:end);
Green   = Green(:,:,1:2:end);
Green = Green(:);
%Red   = Red (:);
outGreen = '\\vault2\Vault2Data\Kelson\B264\18-11-23\OFC_silence\greenchannel.raw';
fwrite(fopen(outGreen, 'w+'),Green,'uint16');
ouRed = '\\vault2\Vault2Data\Kelson\B264\18-11-23\OFC_002\redchannel.raw';
%fwrite(fopen(outGreen, 'w+'),Red,'uint16');