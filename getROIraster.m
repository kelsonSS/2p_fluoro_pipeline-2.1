function raster = getROIraster(spikes,thr,nrois,nf,movtimes,fr)
%spikes -> decon spike struct
%thr -> threshold vals operated on rois
%nrois
%nf -> number of frames
%fr -> frame rate (for lifetime calc)
%nf = numel(find(spikes.chans == 1));
S = zeros(nrois,nf);  %S raster
Sthresh = zeros(nrois,nf);  %S above thresh raster
%loop through cells
for i = 1:nrois    
    %Select cell, i.e., "channel"
    tpos = find(spikes.chans == i);
    %apply threshold and save to Sthresh 
    if ~isempty(tpos)
        %Select spikes that correspond to the selected cell
        vals = spikes.vals(tpos);
        %Select spike values > thr
        tpos2 = find(vals >= thr);
        Sthresh(i,tpos2) = vals(tpos2);        
        S(i,:) = vals;
    end 
end
raster.S = S;
raster.Sthresh = Sthresh;



