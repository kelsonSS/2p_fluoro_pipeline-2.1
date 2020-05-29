





%% init variables 
    Devtotal = zeros(Ncells,Ncells,nmode); DLLtotal = Devtotal;
    wftotal = cell(Ncells,nmode);
    robstotal = cell(Ncells,nmode);
    sig2ftotal = zeros(Ncells,nmode);
    sig2rtotal = zeros(Ncells,Ncells,nmode);
    Bftotal = zeros(Ncells,nmode);
    Brtotal = zeros(Ncells,Ncells,nmode);
    rhatftotal = cell(Ncells,nmode);
    rhatrtotal = cell(Ncells,Ncells,nmode);
    gammaOpttotal = zeros(Ncells,nmode);
    cellidstotal = cell(nmode,1);
 
    %% main loop 
for imd = 1:nmode
    fprintf('Expt %d/%d mode:%d/%d \n', expt,length(dataDir),imd,nmode)
    Resp = respx{imd}(:,:,cellids);
    [N,R,Ncells] = size(Resp);
    %%% GC Analysis: Multi-trial Multi-Cell Mod
    % Estimation Setting:
    LL = 1000;
    LR = 10;
    % Cross-History covariate setting
    WH = 5; %window: average the values over the window over moving window: filters data over 6 Hz that are not calcium signals
    Mhc = 5; %nunber of samples per cell in the model
    
    Whc = [1 WH*ones(1,Mhc-1)]; % Cross-Hist Kernel-make sure that the length latency of interaction + latency of Ca maging
    Lhc = sum(Whc);
    
    % Self-Hist covariate setting
    Mhcs = 5;
    Whcs = [1 WH*ones(1,Mhcs-1)]; % Self-Hist Kernel
    Lhcs = sum(Whcs);
    Mf = (Ncells-1)*Mhc+Mhcs; Mr = (Ncells-2)*Mhc+Mhcs;
    %%% Cross-Validation Step
    gammacv = [0:0.2:3]; % check for 15 diff gamma for each cell
    Whcv = {Whc}; Whscv = {Whcs};
    LLcv = 100; %LL
    LRcv = LR;
    [gammaOpt,Jcost] = CrossValidModTrial2(Resp,gammacv,Whcv,Whscv,LLcv,LRcv);
    %save(['CrossValidData' fileName 'mode' mode{imd} 'Ncells35Wh' num2str(WH) 'Md' num2str(Mhc) 'ModWhc2.mat'],'gammaOpt','gammacv','Whcv','LLcv','LRcv')
    gammaOpttotal(:,imd) = gammaOpt;
    Lm = max([Lhc,Lhcs]);
    Np = N - Lm;
    %%% Form Cross-History Covariate Matrices
    Xcz = zeros(Np,Mhc,Ncells,R);
    for r = 1:R
        Xcz(:,:,:,r) = FormHistMatrixv2(squeeze(Resp(:,r,:)),Whc,Lm); % resp of size [N x Ncells]   
    end
    %%% Form Observation arrays for all cells
    Robs = Resp(Lm+1:end,:,:); %Robs = Resp; 
    % Zero-mean the observation vectors
    % This will remove the need for intercept (mu) estimation
    for cc = 1:Ncells
        Robs(:,:,cc) = Robs(:,:,cc) - ones(Np,1)*mean(Robs(:,:,cc));
    end
    for ct = 1:Ncells % ct: index of target cell
        N_indto = num2str(ct);
        cellC = 1:Ncells; cellC(ct) = [];
        gamma = gammaOpt(ct);
        %%% Form SELF-History Covriate matrix
        Xcsz = zeros(Np,Mhc,R);
       for r = 1:R
            Xcsz(:,:,r) = FormHistMatrixv2(squeeze(Resp(:,r,ct)),Whcs,Lm); % resp of size [N x Ncells]   
        end
        %%% Form Effective response vector
        robsz = Robs(:,:,ct);
        reff = robsz(:); Neff = numel(reff);
        %reffn = reff/sqrt(var(reff));
        %%% Full Model
        %%% Form standardized Full Model Covariate Matrix
        [Xneff,Dnf] = FormFullDesignv3(Xcz,Xcsz,ct);
        %%% Filtering/Estimation Setting
        tic
        [wnf,sig2f,Devf,Bf] = StaticEstim(Xneff,reff,gamma,LL,LR);
     
        %%% Save estimated parameters/arrays
        % Save Tuning AR Coeffs full model for each cell
        wftotal{ct,imd} = wnf;
        sig2ftotal(ct,imd) = sig2f;    
        robstotal{ct,imd} = robsz;
        %Dtotal{ct} = Dnf; % SAVE scalings for cells
        % Reconstructed observation vector full model
        rhatf = zeros(Np,R);
        for r = 1:R
            Xnf = Xneff((r-1)*Np+1:r*Np,:);
            rhatf(:,r) = Xnf*wnf;
        end
        rhatftotal{ct,imd} = rhatf;
        Bftotal(ct,imd) = Bf;
        
        for cf = cellC % cf : index of Electrode/Unit we check causality from to the target 'ct'
            N_indfrom = num2str(cf) ;
            
            
            line_length =  fprintf(['Estimating G-Causality from cell ' N_indfrom ' to cell ' N_indto ' ...']);
            %%% Reduced Model
            Xnefr = Xneff;
            cc = find(cellC == cf);
            Xnefr(:,(cc-1)*Mhc+Mhcs+1:cc*Mhc+Mhcs) = [];
            %%% Estimation Setting           
            [wnr,sig2r,Devr,Br] = StaticEstim(Xnefr,reff,gamma,LL,LR);
        
            % Save estimated params/arrays
            sig2rtotal(ct,cf,imd) = sig2r;
            % Reconstructed observation vector reduced model
            rhatr = zeros(Np,R);
            for r = 1:R
                Xnr = Xnefr((r-1)*Np+1:r*Np,:);
                rhatr(:,r) = Xnr*wnr;
            end
            rhatrtotal{ct,cf,imd} = rhatr;
            Brtotal(ct,cf,imd) = Br;
            %%% Granger Causality metric: Compute Difference of Deviance Statistic
            Devd = Devf - Devr;
            DB = Bf - Br;
            DLL = Devd - DB;
            Devtotal(ct,cf,imd) = Devd;
            DLLtotal(ct,cf,imd) = DLL;
            fprintf(repmat('\b',1,line_length));
        end
    end
    cellidstotal{imd} = cellids;
end

%%% Statistical Test: J-statistics based on FDR
Md = Mhc; %bigger value test higher stat val.
alpha = 0.01; %sig. level; 

DkSmth = Devtotal -Md;
%[GCJL,GCJH] = FDRcontrolBHv4(Devtotal,DkSmth,alpha,Md,wftotal,Whc);
Mlow = 2 %determines which segemetns of Whx are used to dtermine inhibitory vs ecit: limit length of window to avoid phase ambiguoity
[GCJL,GCJH] = FDRcontrolBHv5(Devtotal,DkSmth,alpha,Md,wftotal,Whc,Mlow);

%% Plot Results and Figures

hf = figure;
set(hf, 'position',[10 10 810 610]); 

% find optimal configuration of subplots for most square config 
%( ex ifnmode = 8, rc = [3 3] 
rc = numSubplots(nmode);

for imd = 1:nmode
    cgtxt = {cellidstotal{imd}};
    subplot(rc(1),rc(2),imd), imagesc(Devtotal(:,:,imd))
    cmax = max(max(Devtotal(:,:,imd)));
    cmin = 0;
    
    caxis([cmin cmax])
    title( mode{imd} )
   % title([ mode{imd} ' @ Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd))) ])
    set(gca,'XTick',[1:Ncells],'YTick',[1:Ncells]);
    set(gca,'xticklabel',cgtxt,'yticklabel', cgtxt);%{cellid}
    %xtickangle(90)
    colorbar;
end
%title([ mode{imd} ' @ Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd))) ])
%suptitle(['Estimated deviance @ Md = ' num2str(Md) ])
colormap jet
%saveas(gcf,fullfile(dir_out,['FigDevmap' num2str(Md) 'Whc' num2str(WH) 'gammaAvg' num2str(mean(gammaOpt)) 'Modnnrm.fig']));
saveas(gcf,fullfile(dir_out,['FigDevmap' num2str(Md) 'Whc' num2str(WH)  'Modnnrm.fig']));
close(hf)

cmax = 1; cmin = -1;
hf = figure;
set(hf, 'position',[10 10 810 610]); 
for imd = 1:nmode
    cgtxt = {cellidstotal{imd}};
    subplot(rc(1),rc(2),imd), imagesc(GCJL(:,:,imd))
    caxis([cmin cmax])
    title( mode{imd} )
    %title([ mode{imd} ' @ FDR =' num2str(alpha) ' Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd)))])
    set(gca,'XTick',[1:Ncells],'YTick',[1:Ncells]);
    set(gca,'xticklabel',cgtxt,'yticklabel', cgtxt);%{cellid}
    %xtickangle(90)
    colorbar;
end
%title([ mode{imd} ' @ FDR =' num2str(alpha) ' Whc=' num2str(WH) ' Md= ' num2str(Md) ' gamma= ' num2str(mean(gammaOpttotal(:,imd)))])
%suptitle(['Detected GC maps @ FDR = ' num2str(alpha) ' Md = ' num2str(Md) ])
colormap jet
saveas(gcf,fullfile(dir_out,['FigGCmapsLowLat_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'ModnnrmRun' num2str(irun) '.fig']));
close(hf)



%%% Plot the GC network maps based on Original location of cells
ymax = ceil(max(max(xy))/100)*100 + 10;
%Set min/max for marker size and thickness 
thick_max = 7;
Marker_max = 10;%Marker_min= 10;
aa = 5;
crds = xy;

hf = figure;
set(hf, 'position',[0 0 800 800]); % To set position and scaling of a figure window
for imd = 1:nmode
    XY = xy(cellidstotal{imd},:);
    crdsf = XY; % Set Locations of Selected Cells
    rlen = zeros(Ncells); rtheta = rlen;
    for i = 1:Ncells
        ji = 1:Ncells; ji(i) = [];
        for j = ji
            rlen(i,j) = norm(crdsf(i,:) - crdsf(j,:));
            rtheta(i,j) = atan((crdsf(i,2) - crdsf(j,2))/(crdsf(i,1) - crdsf(j,1)));%abs(crdsf(i,2) - crdsf(j,2))/rlen(i,j);
        end
    end

    hs = subplot(rc(1),rc(2),imd);
    axis square   
    % Plot Markers for All cells first    
    for i = 1:size(crds,1) %v = Channels
        plot(crds(i,1),crds(i,2),'o','MarkerSize',Marker_max,'MarkerFaceColor','w','MarkerEdgeColor','r')
        text(crds(i,1)+5,crds(i,2)+5,num2str(i),'color','r','Fontsize',12)
        hold on
    end
    % Plot Markers for Selec cells    
    pp = Marker_max;
    for i = 1:Ncells %v = Channels
        plot(crdsf(i,1),crdsf(i,2),'o','MarkerSize',pp,'MarkerFaceColor','y','MarkerEdgeColor','r')
        %text(crdsf(i,1)+10,crdsf(i,2)+10,num2str(Channels(i)),'Fontsize',10)
    end
    % Plot GC links
    for r = 1:Ncells
        for c = 1:Ncells
            if c ~= r
                ll= abs(GCJL(r,c,imd))*thick_max ;%+ thick_min
                if ll> thick_max, ll = thick_max; end
                if (ll>0)
                    if r > c, aar = -aa;
                    else aar = aa;
                    end
                    if sign(GCJL(r,c,imd))< 0
                        cll = 'b';
                    else
                        cll = 'r';
                    end
                    offs = aar*[-sin(rtheta(r,c)) cos(rtheta(r,c))];
                    arrow([crdsf(c,1)+offs(1) crdsf(c,2)+offs(2)],[crdsf(r,1)+offs(1) crdsf(r,2)+offs(2)],'Color',cll,'Width',ll,'Length',15,'TipAngle',30);
                end
            end
        end
    end
    axis([0 ymax 0 ymax])
    axis xy
    set(gca,'xticklabel',{[]}, 'yticklabel', {[]});
    ax = gca;
    ax.Color = 'k';
    title(['GC net maps @ mode ' mode{imd}])
end
saveas(gcf,fullfile(dir_out,['FigGCnetworkMaps_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'ModCVModWh2ModnnrmRun' num2str(irun) '.fig']));
close(hf)

    %%% Compute different Statistic measures
    GCnumbers = zeros(nmode,1);
    GClengths = cell(nmode,1);
    GCangles = cell(nmode,1);
    GCavgL = zeros(nmode,1); %GCavgLw = GCavgL;
    
    for imd = 1:nmode
        %%% Compute number of GC links
        GCnumbers(imd) = numel(find(GCJL(:,:,imd)));
        %%% Compute GC link lengths
        GCL = (GCJL(:,:,imd) ~= 0).*rlen;
        GClengths{imd} = GCL(find(GCL));
        GCavgL(imd) = mean(GCL(find(GCL)));
        %GCLw = GCJL(:,:,imd).*rlen;
        %GCavgLw(imd) = mean(GCLw(find(GCLw)));
        %%% Compute GC link angles
        [ii,jj] = find(GCJL(:,:,imd));
        GCtheta = zeros(numel(ii),1);
        for l = 1:numel(ii)
            dy = xy(cellids(ii(l)),2) - xy(cellids(jj(l)),2);
            if dy*rtheta(ii(l),jj(l)) < 0
                GCtheta(l) = rtheta(ii(l),jj(l)) - sign(rtheta(ii(l),jj(l)))*pi;
            else
                GCtheta(l) = rtheta(ii(l),jj(l));
            end
        end
        GCangles{imd} = GCtheta;
    end
    save(fullfile(dir_out,['DataStatsResultsGCanal_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'Modnnrm.mat']),'GCnumbers','GClengths','GCangles','GCavgL','mode') %,'GCavgLw'   

    %%% Plot activation/inhibition effects
    % Plot responses associated with two target cells: Confirm detected GCs visually
    nsubc = 6;
    for imd = 1:nmode
        nGCs = numel(find(GCJL(:,:,imd)))
        nsubr = ceil(nGCs/nsubc); %min(cnt,nGCs)
        
        [Ncells,~,~] = size(Devtotal);
        %cids = cellidstotal{imd};
        Dev = Devtotal(:,:,imd); Devv = Dev(:);
        [~,indev] = sort(Devv,'descend');
        ii = zeros(nGCs,1); jj = ii;
        for l = 1:nGCs
            jj(l) = ceil(indev(l)/Ncells);
            ii(l) = indev(l) - (jj(l)-1)*Ncells;
        end
        GClist = [ii,jj]
        
        ymax = max(max(mean(respx{imd}(:,:,cellids),2)));
        ymin = min(min(mean(respx{imd}(:,:,cellids),2)));        
        hfimd = figure;
        for l = 1:nGCs
            ct = GClist(l,1);
            cf = GClist(l,2);
            cellids = cellidstotal{imd};
            respcf = respx{imd}(:,:,cellids(cf)); respct = respx{imd}(:,:,cellids(ct));
            [N,~] = size(respct);
            subplot(nsubr,min(nsubc,nGCs),l), plot(mean(respct,2))
            hold on
            plot(mean(respcf,2),'r')
            title(['GC(' num2str(cellids(cf)) '\rightarrow' num2str(cellids(ct)) ')'])
            xlim([1,N])
            ylim([ymin ymax])
        end
        saveas(hfimd,fullfile(dir_out,['FigEffectResp_' fileName1 mode{imd} 'mode_Md' num2str(Md) 'Whc' num2str(WH) 'FDR' num2str(alpha) 'ModErlyDlyModnnrmRun' num2str(irun) '.fig']));
        close(hfimd)
    end

    save(fullfile(dir_out,['DataResultsGCanal_' fileName1 '_Md' num2str(Md) 'Whc' num2str(WH) 'ModCVModWh2Modnnrm.mat']),'Devtotal','GCJL','GCJH','wftotal','Whc','Whcs','sig2ftotal','sig2rtotal','robstotal','rhatftotal','rhatrtotal','Bftotal','Brtotal','cellidstotal','gammaOpttotal','LL','LR','alpha')

end 


%% Gather GC data Back into Experiment object 
    Results = [];
    %try
    Results = Collect_GC_Results(dir_out_all);
    %catch
    %   warning('Error in collecting GC results!')
    % end 
%% [Not neccessary to run!] More Detailed Tests on different parts of analysis

%% Check estimated hist kernels of Self + Cross couplings
% imd = 1
% ct = 11
% 
% wf = wftotal{ct,imd};
% 
% figure
% stem(wf)
% 
% numel(find(abs(wf) > std(wf)/5))


%%
% ct = 11
% 
% rct = robstotal{ct,imd};
% 
% figure
% plot(rct)

%% Find Most Significant GC Links
% imd = 1;
% [Ncells,~,~] = size(Devtotal);
% cids = cellidstotal{imd};
% Dev = Devtotal(:,:,imd);
% Devv = Dev(:);
% [~,indev] = sort(Devv,'descend');
% 
% cnt = 20;
% ii = zeros(cnt,1); jj = ii;
% for l = 1:cnt;
%     jj(l) = ceil(indev(l)/Ncells);
%     ii(l) = indev(l) - (jj(l)-1)*Ncells;
% end
% 
% GClist = [ii,jj]
% GClistp = [cids(ii),cids(jj)]

%% Plot responses associated with two target cells: Confirm detected GCs visually
% imd = 1
% ct = 14
% cf = 5
% 
% respcf = robstotal{cf,imd};
% respct = robstotal{ct,imd};
% [Np,~] = size(rhatftotal{1,1});
% 
% figure
% plot(mean(respct,2)/sqrt(var(respct(:))))
% hold on
% plot(mean(respcf,2)/sqrt(var(respct(:))),'r')
% plot(1:Np,mean(rhatftotal{ct,imd},2),'k','linewidth',2)
% plot(1:Np,mean(rhatrtotal{ct,cf,imd},2),'g','linewidth',2)
% title(['Mean activity of cells ' num2str(cf) ' & ' num2str(ct) ' across trials'])
% legend(num2str(ct),num2str(cf))
% xlim([1,Np])
%saveas(gcf,['Sim1SignifGCs_MeanCellAct' num2str(cf) '&' num2str(ct) '.fig']);
%close(gcf)

% iif = 37
% fileName = files(fileIndex(iif)).name
% load(fileName)
%HDFF = Data.eHDFF; % Hit Trial Responses
%MDFF = Data.MDFF; % Miss Trial Responses
%xy = Data.xy; % contain locations of cells
%respx{1} = HDFF; respx{2} = MDFF;
%clear Data HDFF MDFF

% cellids = cellidstotal{imd};
% respcf = respx{imd}(:,:,cellids(cf)); respct = respx{imd}(:,:,cellids(ct));
% [N,~] = size(respct);
% figure
% plot(mean(respct,2))
% hold on
% plot(mean(respcf,2),'r')
% title(['Mean activity of cells ' num2str(cf) ' & ' num2str(ct) ' across trials'])
% legend(num2str(ct),num2str(cf))
% xlim([1,N])
% 
% figure
% plot(respct,'b')
% hold on
% plot(respcf,'r')
% title(['Activity of cells ' num2str(cf) ' & ' num2str(ct) ' within all trials'])
% xlim([1,N])
% saveas(gcf,['Sim1SignifGCs_CellAct' num2str(cf) '&' num2str(ct) '.fig']);
%close(gcf)

%%
% imd = 2
% ct = 10
% cf = 11
% 
% rhatct = rhatftotal{ct,imd};
% rhatcf = rhatftotal{cf,imd};
% 
% figure
% subplot(1,2,1),plot(rhatct,'b'), hold on, plot(mean(rhatct,2),'k','linewidth',2)
% title(['Activity across trials @' mode{imd} ' ct=' num2str(cellids(ct))])
% xlim([0 N-Lm])
% subplot(1,2,2),plot(rhatcf,'r'), hold on, plot(mean(rhatcf,2),'k','linewidth',2)
% title(['Activity across trials @' mode{imd} ' cf=' num2str(cellids(cf))])
% xlim([0 N-Lm])
% 
% rr = 9
% figure
% subplot(1,2,1),plot(rhatct(:,rr),'b'), hold on, plot(mean(rhatct,2),'k','linewidth',2)
% title(['Activity @' mode{imd} ' ct=' num2str(cellids(ct))])
% xlim([0 N-Lm])
% subplot(1,2,2),plot(rhatcf(:,rr),'r'), hold on, plot(mean(rhatcf,2),'k','linewidth',2)
% title(['Activity @' mode{imd} ' cf=' num2str(cellids(cf))])
% xlim([0 N-Lm])
% 
% %%
% ct = 16
% wnf = wftotal{ct,imd};
% 
% figure
% stem(wnf)
% 
% %% Plot the sample figure of activations 
% hf = figure;
% l = [1,2,4,8]
% for i = 1:4
%     ct = GClist(l(i),1);
%     cf = GClist(l(i),2);
%     cellids = cellidstotal{imd};
%     respcf = respx{imd}(:,:,cellids(cf)); respct = respx{imd}(:,:,cellids(ct));
%     [N,~] = size(respct);
%     subplot(1,4,i), plot(mean(respct,2))
%     hold on
%     plot(mean(respcf,2),'r')
%     title(['GC(' num2str(cellids(cf)) '\rightarrow' num2str(cellids(ct)) ')'])
%     xlim([1,N])
%     ylim([-0.7 +0.7])
% end
% 



%% PreProcessing: Check response traces
%% Plot Hit/Miss Trial Activties of specific cell
% cc = 1;
% imd = 1;
% respcc = squeeze(resp{imd}(:,:,cc));
% 
% figure
% plot(respcc(:,:))
% hold on
% plot(mean(respcc,2),'r','linewidth',2)
% axis tight
% title(['Response during Hits for All Trials cell #' num2str(cc)])

%%% Plot Hit/Miss Activties of all cells within trial
% rr = 15;
% imd = 1;
% resprr = squeeze(respx{imd}(:,rr,:));
% %resprr = squeeze(resp{imd}(:,rr,cellids));
% 
% figure
% plot(resprr(:,:))
% axis tight
% title(['Response during Hits for All cells within trial #' num2str(rr)])

% %%
% Nc = 6;
% figure
% for i = 1:Nc
%     subplot(Nc,1,i), plot(respcc(:,i))
%     axis tight
% end