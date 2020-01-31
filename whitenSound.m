sf = 200000; % sampling frequency

W_time = length(WhiteningSpec) / sf ; %length of whiteningSpec in seconds 
stim_time = length(stim)/ sf; % length of stimulus in seconds 

chunk_num_iter = ceil(stim_time / W_time);
chunk_size = length(WhiteningSpec);
stimWhite = zeros(1,length(stim));
for chunk = 1:chunk_num_iter
    fprintf('chunk: %d of %d \n',chunk, chunk_num_iter)
    cp = [1 chunk_size] + chunk_size * (chunk-1); % chunk pointers
    
    if chunk == chunk_num_iter
        cp(2) = length(stim);
       chunk_size = cp(2) - cp(1) +1;
    end 
    

    stimSpec = fft(stim(cp(1):cp(2)),length(WhiteningSpec));
    
stimSpecPhase = angle(stimSpec);
stimSpecdB=VolumeConversion(abs(stimSpec),'V2dB',mic);
WhiteningSpecdB=VolumeConversion(WhiteningSpec,'V2dB',mic);
stimSpecdBWhite = stimSpecdB + WhiteningSpecdB;
stimSpecWhite = VolumeConversion(stimSpecdBWhite,'dB2V',mic);
stimSpecWhite = stimSpecWhite.*exp(j.*stimSpecPhase);
stimWhite_temp  = real(ifft(stimSpecWhite));
stimWhite(cp(1):cp(2)) = stimWhite_temp(1:chunk_size); 

end 
audiowrite('\\vault3\data\kelson\WhiteNoise.wav',stimWhite,sf)