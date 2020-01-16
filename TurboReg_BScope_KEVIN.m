addpath('/applications/fiji.app/scripts')
Miji;

root = '/Volumes/Vault2Data/Kevin/Tone_Rearing/PV-cre/';
% folder_name = {'11082016_m810_longtone_lowA1','11082016_m8130_longtone_lowAAF',...
%     '11082016_m8131_longtone_lowA2'};
% folder_name = { 'stimAq08001','stimAq09','stimAq10','stimAq11','stimAq12',...
%     'stimAq13','stimAq14','stimAq15'};
% folder_name = {'11262016_m721_longtone_AAF',...
%    '11262016_m723_longtone_highA1','11262016_m962_longtone_AAF'};
  folder_name = {     
        'Harvey/layerIV/unknown',
        'Ben/layerIV/unknown'
     };
     ;

N = length( folder_name);
target_img_name = repmat({'AVG_greenchannelregistered.tif'},N,1);
% target_img_name = {'AVG_redChan.tif','AVG_redChan.tif',...
%     'AVG_redChan.tif','AVG_redChan.tif'};
nTotal = 100000 * ones(N,1);
width = 512*ones(N,1);
height = 512*ones(N,1);
% dual_chan = [2,1;2,1;2,1;2,1];
%dual_chan = repmat([2,1],N,1); % put back if red channel exists
dual_chan = [nan,nan];
correct_which_channel = repmat({'green'},N,1);
% if single channel then put dual_chan = [nan,nan]; if dual channel recording,
% in each row, the first number indicates which channel is green, second number
% indicates which is red
useRed2CorrectGreen = 0;

%% cycle through all folders
nFramePerLoop = 5000; % # of images to export in each loop

for i = 1:length(folder_name)%1:2
    tic
    fldn = folder_name{i};
    
    fprintf('Processing folder %s\n',fullfile(root,fldn));
    rawfile = dir(fullfile(root,fldn,'Image_0001_0001.raw'));
    %rawfile = dir(fullfile(root,fldn,'greenchannel.raw'));
    
    % the N by 2 matrix to store shifts
    shift = zeros(nTotal(i),2);
    
    % path to save corrected images
    path2save = fullfile(root,fldn,'shifted_img');
    if ~isdir(path2save)
        mkdir(path2save);
    end
    % the path to save original images
    if ~useRed2CorrectGreen
        ori_img_sfldn = [correct_which_channel{i} 'Chan'];
        path2save_ori = fullfile(root,fldn,ori_img_sfldn);
        if ~isdir(path2save_ori)
            mkdir(path2save_ori);
        end
    else
        path2save_ori = fullfile(root,fldn,'redChan');
        if ~isdir(path2save_ori)
            mkdir(path2save_ori);
        end
        path2save_ori = fullfile(root,fldn,'greenChan');
        if ~isdir(path2save_ori)
            mkdir(path2save_ori);
        end
    end
    % the gap to read images depend on if it is dual channel
    if isnan(dual_chan(i,1)) % single channel
        gap = 0;
    else % dual channel
        gap = 3 * width(i)*height(i)*2;
    end
    
    
    last = 0;rng = [];
    while last < nTotal(i)
        tmp = min(nFramePerLoop,nTotal(i)-last);
        rng = [rng; last + 1, last + tmp];
        last = last + tmp;
    end
    nTotalImage = 0;
    
    switch correct_which_channel{i}
        case 'green'
            which_channel = 1;
        case 'red'
            which_channel = 2;
    end
    
    % generate .mat file to save the shift information for each frame
    shift_x = zeros(nTotal(i),1);
    shift_y = zeros(nTotal(i),1);
    landmark_path = fullfile(root,fldn,'landmarks.mat');
    save(landmark_path,'shift_x','shift_y','-v7.3');  % preallocate shift info
    hMat = matfile(landmark_path,'Writable',true);
    
    for k = 1:size(rng,1)
        % offset also depend on dual or single channel
        if isnan(dual_chan(i,1)) % single channel
            offset = 0 + (k-1)*nFramePerLoop*width(i)*height(i)*2; % offset in bytes, 16bit = 2bytes
        else
            offset = (dual_chan(i,which_channel)-1)*width(i)*height(i)*2  + (k-1)*nFramePerLoop*width(i)*height(i)*2*4; % offset in bytes, 16bit = 2bytes
        end
        if offset >= rawfile.bytes
            break;
        end

        n_img_2_load = min( ( rawfile.bytes / (gap + width(i)*height(i)*2) - (k-1)*nFramePerLoop )  , nFramePerLoop );
        nTotalImage = nTotalImage + n_img_2_load;
        
        fprintf('Processing image %d to %d\n',rng(k,1),rng(k,1) + n_img_2_load - 1);

        % MUST load target image first
        MIJ.run('Open...',sprintf('open=%s',fullfile(root,fldn,target_img_name{i})));
        
        % load images from .raw file
        MIJ.run('Raw...',...
            sprintf('open=%s image=[16-bit Unsigned] width=%d height=%d offset=%d number=%d gap=%d little-endian',...
            fullfile(root,fldn,'Image_0001_0001.raw'),width(i),height(i),offset,diff(rng(k,:))+1,gap ));
        
        % can export tiff files here actually
        if ~useRed2CorrectGreen
            MIJ.run('Image Sequence... ', sprintf('format=TIFF name=%s start=%d digits=5 save=%s ',...
                correct_which_channel{i},rng(k,1),fullfile(path2save_ori,sprintf('%s%05d.tif',correct_which_channel{i},rng(k,1)))));
        else
            MIJ.run('Image Sequence... ', sprintf('format=TIFF name=%s start=%d digits=5 save=%s ',...
                'red',rng(k,1),fullfile(root,fldn,'redChan',sprintf('%s%05d.tif','red',rng(k,1)))));
            
            % load green channel
            offset_green = (dual_chan(i,1)-1)*width(i)*height(i)*2  + (k-1)*nFramePerLoop*width(i)*height(i)*2*4;
            MIJ.run('Raw...',...
                sprintf('open=%s image=[16-bit Unsigned] width=%d height=%d offset=%d number=%d gap=%d little-endian',...
                fullfile(root,fldn,'Image_0001_0001.raw'),width(i),height(i),offset_green,diff(rng(k,:))+1,gap ));
            MIJ.run('Image Sequence... ', sprintf('format=TIFF name=%s start=%d digits=5 save=%s ',...
                'green',rng(k,1),fullfile(root,fldn,'greenChan',sprintf('%s%05d.tif','green',rng(k,1)))));
            MIJ.run('Close');
            
        end
                
        % call TurboReg
        path_2_save_landmarks = fullfile(root,fldn);
        if ~ismac
            myTurboRegOption = sprintf('-RunBatch %s\\', path_2_save_landmarks);
        else
            myTurboRegOption = sprintf('-RunBatch %s/', path_2_save_landmarks);
        end
        tic
        MIJ.run('TurboReg ', myTurboRegOption);
        
        % need to wait before TurboReg finishes correcting this folder
        %     we do this by constantly checking whether a 'finished.txt' exists in
        %     root\folder_name\
        while ~exist(fullfile(root,fldn,'finished.txt'),'file')
            pause(10);
        end

        delete(fullfile(root,fldn,'finished.txt')); % after we detect it is finished, delete the txt file for next loop 

        t = toc;
        fprintf('FINISHED correcting %d images in %d minute(s)\n',n_img_2_load,round(t/60));
        MIJ.run('Close All');
       
        % load landmarks.txt and write info into landmarks.mat
        
        landmark = zeros(n_img_2_load,4);
        
        fid = fopen(fullfile(root,fldn,'landmarks.txt'),'r');
        Transformation = fscanf(fid,'%s\n',1);
        src_sz = fscanf(fid,'\nSource size\n%d %d\n',[1,2]);
        tgt_sz = fscanf(fid,'\nTarget size\n%d %d\n',[1,2]);
        j=1;
        while ~feof(fid)
            landmark(j,:) = fscanf(fid,'\nRefined source landmarks\n%f %f\n\nTarget landmarks\n%f %f\n',[1,4]);
            j = j+1;
        end
        shift_x_tmp = landmark(:,3) - landmark(:,1);
        shift_y_tmp = landmark(:,4) - landmark(:,2);
        hMat.shift_x(rng(k,1):rng(k,1) + n_img_2_load -1 ,1) = shift_x_tmp; % save to the matfile
        hMat.shift_y(rng(k,1):rng(k,1) + n_img_2_load -1,1) = shift_y_tmp;
        
        shift_x(rng(k,1):rng(k,1) + n_img_2_load -1,1) = shift_x_tmp; % also save to workspace
        shift_y(rng(k,1):rng(k,1) + n_img_2_load -1,1) = shift_y_tmp;
        fclose(fid);
        
    end
    % finish correcting all images
    % Do actual registration and write shifted img
    path_2_save = fullfile(root,fldn,'shifted_img');
    if ~isdir(path_2_save)
        mkdir(path_2_save);
    end
    shift_x = shift_x(1:nTotalImage,:); % truncate at the # of total images
    shift_y = shift_y(1:nTotalImage,:);
    
    hMat.shift_x = shift_x;
    hMat.shift_y = shift_y;
    if ~useRed2CorrectGreen
        sample = imread(fullfile(root,fldn,ori_img_sfldn,sprintf('%s%05d.tif',correct_which_channel{i},1)));
    else
        sample = imread(fullfile(root,fldn,'greenChan',sprintf('%s%05d.tif','green',1)));
    end
    avg_img = zeros(size(sample));
    
    fprintf('Saving images...\n')
    tic
    parfor k = 1:nTotalImage
        if useRed2CorrectGreen
           % always use green channel in this case
           old_im = im2double(imread(fullfile(root,fldn,'greenChan',sprintf('green%05d.tif',k))));
        else
           old_im = im2double(imread(fullfile(root,fldn,ori_img_sfldn,sprintf('%s%05d.tif',correct_which_channel{i},k))));    
        end
        new_im = imtranslate(old_im,[shift_x(k),shift_y(k)],'cubic','OutputView','same');
        avg_img = avg_img+new_im;
        imwrite(im2uint16(new_im),fullfile(path_2_save,sprintf('Registered%05d.tif',k)));
    end
    avg_img = avg_img / nTotalImage;
    imwrite(im2uint16(avg_img),fullfile(root,folder_name{i},'AVG_shifted_img_TurboReg.tif'));
    t = toc;
    fprintf('FINISHED generating corrected images, time spent: %d minutes\n',round(t/60));

end





