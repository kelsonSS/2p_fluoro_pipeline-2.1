function num_frames =  FindRawImgSize(fileID,img_dimensions,byte_size)
% this function takes a file handle that has been opened by fopen 
% along with the size of how many bytes per number and expected single image size
% and returns number of expected frames in the movie
%
% 
% uint16 = 2 bytes 
% Kelson shilling-scrivo 2020




if ~exist('byte_size','var')
    byte_size = 2;
end 

if ~exist('img_dimensions','var')
    img_dimensions = [512 512];
end 


% move pointer to end of file
fseek(fileID,0,'eof');

          % number of images = (# of bytes/ bytes per number)  /  (numbers per image) 
num_frames = (ftell(fileID)/byte_size) / (img_dimensions(1) *img_dimensions(2));

% rewind pointer to begging of frame
fseek(fileID,0,'bof');



end 

