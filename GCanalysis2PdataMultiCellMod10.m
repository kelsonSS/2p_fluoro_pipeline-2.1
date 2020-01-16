%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Patrik Data GC Analysis
%%%% Oct 19th 2016 : PairWise Multi-Trial GC Analysis
%%%% Oct 25th 2016 : MultiCell GC Analysis
%%%% Nov 2nd 16: Multi-resolution History dependence
%%%% Nov 4th 16: The integrative AR models (didn't work well)
%%%% Dec 1st 16: Change of Standardization/Normalization procedure on
%%%% repetitions (Trial-dependent variance scenario)
%%%% Dec 2nd 16: Modify the new resolved case to Trial-independent variance scenario
%%%% Dec 23rd 16: Sim10 modification of code based on new iterations for
%%%% Uniform Sparsity Level 


cd('C:\Users\Kelson\Dropbox\GCanalysis2PdataMultiCellMod10')


load('m509_4secTone.mat')
load('stim_seq_m509_4secTone.mat')
load('mask_m509_4secTone.mat')
load('m509_4secTone_param.mat')

%% Preprocessing Cell activity
resp = traces_cell - 0.7*traces_neuropil;

% Prepare stim sequence of onset-offsets
Fs = 30; % Sampling Freq (Hz)
Tonset = 4; % Duration of onset tones (sec)
Nonset = Tonset*Fs; % Number of Onset samples
framesOnset = find(freq_seq); 
Ns = framesOnset(2)-framesOnset(1); % # of samples during each Onset-Offset

StimSeq = zeros(size(freq_seq));
for i = 1:numel(framesOnset)
    StimSeq(framesOnset(i):framesOnset(i)+Nonset-1) = freq_seq(framesOnset(i))*ones(1,Nonset);
end

% Choose values: [1:9] for 9 different frequencie tones 
ToneSeq = zeros(size(freq_seq));
ToneSeq(framesOnset) = 2*log2(freq_seq(framesOnset)/4)+1; 
Ntone = max(ToneSeq);

%% Filtering & Estimation
%%% Multi-trial Analysis

%%% Cell Selection: Mod 1
[vr,inds] = sort(var(resp'),'descend');
Ncells = 10;
cellid = inds(1:Ncells);
cellid = sort(cellid);

% Complete the last trial of last tone: zero-pad
%Resp = [resp(cellid,:), zeros(Ncells,max(find(ToneSeq > 0))+Ns-1 - size(resp,2))];
% Remove the last incomplete trial of last tone
Resp = resp(cellid,1:max(find(ToneSeq > 0))-1); ToneSeq(max(find(ToneSeq > 0))) = 0;

%%% GC Anlysis ONSET
% Choose 'ON' for Onset and 'OFF' for Offset
mode = 'OFF' 

Devtotal = zeros(Ntone,Ncells,Ncells);
wtotal = cell(Ntone,Ncells);
Dtotal = cell(Ntone,Ncells);
SNR = zeros(Ntone,Ncells);
rhatf = cell(Ntone,Ncells);
robst = cell(Ntone,Ncells);
rhatr = cell(Ntone,Ncells,Ncells);
sig2rt = cell(Ntone,Ncells,Ncells);
sig2ft = cell(Ntone,Ncells);
Bftotal = zeros(Ntone,Ncells);
Brtotal = zeros(Ntone,Ncells,Ncells);

%%% Prepare Measurement Data Matrices X(k)
% Choose order of AR filters for sparse recovery 

%%% Cross-history kernel (multi-resolution window)
WhOpt = [1 1 1 2 2 3]; 
Whc = WhOpt; Mhc = numel(Whc);
Lhc = sum(Whc);
%%% Self-history kernel
Whcs = WhOpt; Mhcs = numel(Whcs);
Lhcs = sum(Whcs);
Lm = max([Lhc,Lhcs]);

% Filtering Setting:
LL = 1000; LR = 20;
gamma = 1;

for Tid = 1:Ntone 
    disp(['Computing GC links for Tone #' num2str(Tid)])
    if strcmp(mode,'ON')
        N = Nonset; 
        indseq = find(ToneSeq == Tid); R = numel(indseq);
    elseif strcmp(mode,'OFF')
        N = Ns-Nonset;
        indseq = find(ToneSeq == Tid)+Nonset; R = numel(indseq);
    else
        disp('WRONG INPUT MODE...')
        break
    end
    
    %%% Extract the cell response trials to the tone Tid
    Xtone = zeros(N+Lm,R,Ncells);
    for i = 1:Ncells
        for r = 1:R
            Xtone(:,r,i) = Resp(i,indseq(r)-Lm:indseq(r)+N-1)'; 
        end
    end
    
    %%% Form the cross-history covariates based on cross-hist kernel
    X = zeros(N,Mhc,R,Ncells);
    for i = 1:Ncells
        for r = 1:R
            for k = 1:N
                for j = 1:Mhc
                    X(k,Mhc-j+1,r,i) = mean(Xtone(k+sum(Whc(Mhc-j+2:Mhc)):k+sum(Whc(Mhc-j+1:Mhc))-1,r,i));
                end
            end
        end
    end
    % Zero-mean the cross-history covariates before regression and filtering:
    % [NECESSARY in terms of (X'X)/N approx= I orthonormal conditions]
    for i = 1:Ncells
        for r = 1:R
            X(:,:,r,i) = X(:,:,r,i) - ones(N,1)*mean(X(:,:,r,i));
        end
    end
    
    %%% Form Observation vectors for all cells
    Robs = Xtone(Lm+1:end,:,:);
    % Zero-mean the observation vectors
    for i = 1:Ncells
        Robs(:,:,i) = Robs(:,:,i) - ones(N,1)*mean(Robs(:,:,i));
    end
    
    Mf = (Ncells-1)*Mhc+Mhcs; % number of parameters full model
    Mr = (Ncells-2)*Mhc+Mhcs; % number of parameters reduced model
    
    for ct = 1:Ncells % ct: index of target cell
        N_indto = num2str(cellid(ct));
        robsz = Robs(:,:,ct);
        cellC = 1:Ncells; cellC(ct) = [];
        
        %%% Form SELF-History Covriates
        Xcs = zeros(N,Mhcs,R);
        for r = 1:R
            for k = 1:N
                for j = 1:Mhcs
                    Xcs(k,Mhcs-j+1,r) = mean( Xtone(k+sum(Whcs(Mhcs-j+2:Mhcs)):k+sum(Whcs(Mhcs-j+1:Mhcs))-1,r,ct) );
                end
            end
        end
        % Zero-mean Self-History Components
        for r = 1:R
            Xcs(:,:,r) = Xcs(:,:,r) - ones(N,1)*mean(Xcs(:,:,r));
        end
        
        %%% FULL Model
        % Form Full model design/covariate matrix Xf
        Xfc = Xcs;
        for cc = 1:Ncells-1
            Xfc = cat(2,Xfc,X(:,:,:,cellC(cc)));
        end
        Xf = Xfc; 

        tic
        %%% Filtering/Estimation Setting
        % Set initials
        w0f = zeros(Mf,1); wl = w0f;
        %%% FORM an EFFECTIVE covariate matrix
        % Standardize the effective covariate matrix Xeff
        Xeff = zeros(R*N,Mf); 
        for r = 1:R
            Xeff((r-1)*N+1:r*N,:) = Xf(:,:,r);
        end
        Dnf = diag(sqrt(var(Xeff)));
        Xneff = Xeff/Dnf;
        % Effective response vector
        reff = robsz(:); Neff = numel(reff);
        
        %%% NEW ITERATIVE UPDATE PROCEDURE 
        XXneff = Xneff'*Xneff; % Effective autocovariance matrix
        Xdneff = Xneff'*reff; % Effective cross-covariance vector
        % Initial estimate of the variance: ML estimate
        wML = XXneff\Xdneff;
        sig2f0 = (norm(reff - Xneff*wML)^2)/Neff; sig2fl = sig2f0;
        for l = 1:LR
            XXfl = XXneff/sqrt(sig2fl*Neff);
            Xdfl = Xdneff/sqrt(sig2fl*Neff);
            % Step size selection
            rhom = max(abs(eig(XXfl))); al = 0.9/rhom;
            % Perform Proximal gradiant iterations
            for ll = 1:LL                
                gfl = (Xdfl - XXfl*wl); % Form gradient 
                wl = SoftThreshold(wl + al*gfl ,gamma*al); 
            end
            % Update variance based on the updated wl estimate
            sig2fl = (norm(reff - Xneff*wl)^2)/Neff; % Trial-independent version
        end
        wnf = wl; gf = gfl;        
        sig2f = sig2fl;
        % Compute Unbiased Deviance
        sig2ft{Tid,ct} = sig2fl;
        Bf = gf'*(XXfl\gf)*sqrt(Neff/sig2f);
        Devf = -Neff*log(sig2f) + Bf;
        toc
        
        %%% Save estimated parameters/arrays
        % Save Tuning AR Coeffs full model for each cell
        wtotal{Tid,ct} = wnf;
        % SAVE scalings for different Trials
        Dtotal{Tid,ct} = Dnf;
        % Reconstructed observation vector full model
        rhati = zeros(N,R);
        Xnf = zeros(size(Xf));
        for r = 1:R
            Xnf(:,:,r) = Xf(:,:,r)/Dnf;
            rhati(:,r) = Xnf(:,:,r)*wnf;
        end
        rhatf{Tid,ct} = rhati;
        robst{Tid,ct} = robsz;
        Bftotal(Tid,ct) = Bf;
        
        for cf = cellC % cf : index of Electrode/Unit we check causality from to the target 'ct'
            N_indfrom = num2str(cellid(cf)) ;
            display(['Estimating Causality from cell ' N_indfrom ' to cell ' N_indto ' ...'])

            %%% Reduced Model
            Xr = Xf;
            cc = find(cellC == cf);
            Xr(:,(cc-1)*Mhc+Mhcs+1:cc*Mhc+Mhcs,:) = [];

            %%% Estimation Setting           
            w0r = zeros(Mr,1); wl = w0r;
            % Standardize the effective covariate matrix Xeff
            Xefr = zeros(R*N,Mr); 
            for r = 1:R
                Xefr((r-1)*N+1:r*N,:) = Xr(:,:,r);
            end
            Dnr = diag(sqrt(var(Xefr)));
            Xnefr = Xefr/Dnr;            
            %%% NEW ITERATIVE UPDATE PROCEDURE
            XXnefr = Xnefr'*Xnefr;
            Xdnefr = Xnefr'*reff;
            % Initial estimate of the variance: ML estimate
            wML = XXnefr\Xdnefr;
            sig2r0 = (norm(reff - Xnefr*wML)^2)/Neff;
            sig2rl = sig2r0;
            for l = 1:LR
                XXrl = XXnefr/sqrt(sig2rl*Neff);
                Xdrl = Xdnefr/sqrt(sig2rl*Neff);
                % Step size selection
                rhom = max(abs(eig(XXrl))); al = 0.9/rhom;
                % Perform Proximal gradiant iterations
                for ll = 1:LL
                    grl = (Xdrl - XXrl*wl); % Form gradient 
                    wl = SoftThreshold(wl + al*grl ,gamma*al); 
                end
                % Update variance based on updated wl estimate
                sig2rl = (norm(reff - Xnefr*wl)^2)/Neff; % Trial-independent version
            end
            wnr = wl; gr = grl;        
            sig2r = sig2rl;
            Br = grl'*(XXrl\grl)*sqrt(Neff/sig2r);
            Devr = -Neff*log(sig2r) + Br ;
            % Save estimated params/arrays
            sig2rt{Tid,ct,cf} = sig2r;
            % Reconstructed observation vector reduced model
            rhatj = zeros(N,R);
            Xnr = zeros(size(Xr));
            for r = 1:R
                Xnr(:,:,r) = Xr(:,:,r)/Dnr;
                rhatj(:,r) = Xnr(:,:,r)*wnr;
            end
            rhatr{Tid,ct,cf} = rhatj;
            Brtotal(Tid,ct,cf) = Br;
            
            %%% Granger Causality metric: Compute Difference of Deviance Statistic
            Dn = (Devf - Devr);
            Devtotal(Tid,ct,cf) = Dn; 
        end
    end
end

%% Statistical Test: FDR Control based on BH
alpha = 0.01 % Significance level
Nl = Ncells*(Ncells-1); % Number of links
pth = alpha*(1:Nl)/(Nl*log(Nl));
alpham = alpha*(Nl+1)/(2*Nl*log(Nl));
Md = Mf - Mr % difference dimensionality

GCJ = zeros(size(Devtotal));
% Compute p-values
Pv = chi2cdf(Devtotal,Md,'upper');

for t = 1:size(Pv,1)
    Pvt = squeeze(Pv(t,:,:));
    [Pvtsort, Indsrt] = sort(Pvt(:));
    Pvtsortx = Pvtsort(1:Nl);
    cnt = 1;
    while Pvtsortx(cnt) < pth(cnt)
        cnt = cnt+1;
    end
    siglind = Indsrt(1:cnt-1);
    cols = ceil(siglind/Ncells);
    rows = siglind - (cols-1)*Ncells;
    for ee = 1:numel(siglind)
        % Compute J-statistics using raw estimate of non-centrality equal to nu = Dev - Md
        GCJ(t,rows(ee),cols(ee)) = 1 - alpham - ncx2cdf( chi2inv(1-alpham,Md),Md, Devtotal(t,rows(ee),cols(ee)) -Md ); % 1 - Pv(t,rows(ee),cols(ee));
    end
end

%%% Add Effective Sign Excit/Inhib Nature
GCJs = zeros(size(GCJ));
for Tid = 1:Ntone
    GCt = squeeze(GCJ(Tid,:,:));
    [ii,jj] = find(GCt);
    GCsgn = zeros(size(GCt));
    for ll = 1:numel(ii)
        whati = wtotal{Tid,ii(ll)}; 
        CellC = 1:Ncells; CellC(ii(ll))=[];
        jji = find(CellC == jj(ll));
        GCsgn(ii(ll),jj(ll)) = sign(Whc*whati((jji-1)*Mhc+Mhcs+1:jji*Mhc+Mhcs));
    end
    GCJs(Tid,:,:) = GCsgn.*GCt;
end

%%%
cmax = max(max(max(Devtotal)));
cmin = 0;
figure
for Tid = 1:Ntone
    subplot(3,3,Tid), imagesc(squeeze(Devtotal(Tid,:,:)))
    caxis([cmin cmax])
    title(['f = ' num2str(2.^((Tid-1)/2)*4) 'kHz'])
    set(gca,'XTick',[1:Ncells],'YTick',[1:Ncells]);
    set(gca,'xticklabel',{cellid},'yticklabel', {cellid});
    %xtickangle(90)
end
colormap jet
suptitle(['Deviance estimates @ ' mode ' mode, Nhc= [' num2str(Whc) '], N= ' num2str(N)])

cmax = 1; cmin = -1;
figure
for Tid = 1:Ntone
    subplot(3,3,Tid), imagesc(squeeze(GCJs(Tid,:,:)))
    caxis([cmin cmax])
    title(['f = ' num2str(2.^((Tid-1)/2)*4) 'kHz'])
    set(gca,'XTick',[1:Ncells],'YTick',[1:Ncells]);
    set(gca,'xticklabel',{cellid},'yticklabel', {cellid});
    %xtickangle(90)
end
colormap jet
suptitle(['J-statistics GC @ ' mode ' mode, Nhc= [' num2str(Whc) '], FDR= ' num2str(alpha) ', N= ' num2str(N)])

%% PLOT the reconstructed responses using full & reduced models
% Choose different Tones and GC links (cf ==> ct)
Tid = 5
ct = 1
cf = 5

R = numel(find(ToneSeq == Tid));
figure
for  r = 1:R
    subplot(5,2,r)
    plot(robst{Tid,ct}(:,r))
    hold on
    plot(rhatf{Tid,ct}(:,r),'r')
    plot(rhatr{Tid,ct,cf}(:,r),'g')
    xlim([1 N])
end
suptitle(['Reconstructed responses for all trials for Tid=' num2str(Tid) ' , ct=' num2str(ct) ' , cf=' num2str(cf) ' , Md=' num2str(Md) ' , gamma=' num2str(gamma)])

figure
for  r = 1:R
    subplot(5,2,r)
    plot(robst{Tid,ct}(:,r))
    hold on
    plot(robst{Tid,cf}(:,r),'r')
    xlim([1 N])
end
suptitle(['Trial Responses for Tid=' num2str(Tid) ' , ct=' num2str(ct) ' , cf=' num2str(cf) ' , Md=' num2str(Md) ' , gamma=' num2str(gamma)])
%saveas(gcf,['Sim8GCmapsCellSel1_' mode 'set_Rev9CVMd' num2str(Md) 'gamma'  num2str(gamma) '_RecTrialsTid' num2str(Tid) 'GC' num2str(cf) 'to' num2str(ct) '.fig']);

Devtotal(Tid,ct,cf)

%% PLOT difference of full & reduced models for all trials
% Choose different Tones and GC links (cf ==> ct)
Tid = 2
ct = 9
cf = 2

R = numel(find(ToneSeq == Tid));
figure
for  r = 1:R
    subplot(5,2,r)
    plot(rhatf{Tid,ct}(:,r)-rhatr{Tid,ct,cf}(:,r),'r')
    xlim([1 N])
    ylim([-0.005 +0.005])
end
suptitle(['Reconstructed Trials for Tid=' num2str(Tid) ' , ct=' num2str(ct) ' , cf=' num2str(cf) ' , Md=' num2str(Md) ' , gamma=' num2str(gamma)])


%% Cross-History Coefficients Analysis
%%% Explore the estimated cross-history coeffs
% Specify the Tone number Tid

Tid = 5
GCt = squeeze(GCJs(Tid,:,:));
[ii,jj] = find(GCt);
nlnk = numel(ii);
if nlnk > 0
    nsubc = 4; nsubr = ceil(nlnk/nsubc); 
    figure
    for ll = 1:nlnk
        ct = ii(ll); cf = jj(ll);
        wtnct = wtotal{Tid,ct};
        % Find Coeffs associated with cell cf
        cellC = 1:Ncells; cellC(ct) = [];
        cc = find(cellC == cf);
        wtnctcf = wtnct((cc-1)*Mhc+Mhcs+2:cc*Mhc+Mhcs+1);

        subplot(nsubr,nsubc,ll), stem(wtnctcf)
        title(['ct=' num2str(ct) ' , cf=' num2str(cf)])
    end
    suptitle(['Estimated coeffs of significant GC links at Tone f=' num2str(2.^((Tid-1)/2)*4)])
else
    disp(['No GC links for tone #' num2str(Tid)])
end

%% PLOT NETWORK Causal MAPs

% Set min/max for marker size and thickness 
thick_max = 5;
Marker_max = 10;
ymax = 510 
crdsf = center(cellid,:);
rlen = zeros(Ncells);
rtheta = rlen;
for i = 1:Ncells
    ji = 1:Ncells; ji(i) = [];
    for j = ji
        rlen(i,j) = norm(crdsf(i,:) - crdsf(j,:));
        rtheta(i,j) = abs(crdsf(i,2) - crdsf(j,2))/rlen(i,j);
    end
end

drawArrow = @(x,y,lw,c)quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0,'LineWidth',lw,'color',c);
aa = 7;

h1 = openfig('cells_m509_4secTone.fig','reuse');% open figure
ax1 = gca; % get handle to axes of figure
fig1 = get(ax1,'children'); %get handle to all the children in the figure

for Tid = 1:Ntone

    hf = figure;
    set(hf, 'position',[0 0 1000 1000]); % To set position and scaling of a figure window

    s1 = subplot(10,1,[1 10]);
    ax = get(s1,'Position');
    set(s1,'Position',ax);
    copyobj(fig1,s1); % creates copies of graphics objects and their descendants and assigns the objects to the new parent.
    axis square 
    hold on
    colormap gray

    % Plot markers     
    pp = Marker_max;
    for i = 1:Ncells 
        plot(crdsf(i,1),crdsf(i,2),'o','MarkerSize',pp,'MarkerFaceColor','y','MarkerEdgeColor','r')
    end

    % Plot connections
    for r = 1:Ncells
        for c = 1:Ncells
            if c ~=r
                ll = abs(GCJs(Tid,r,c))*thick_max ;%+ thick_min
                if ll> thick_max 
                    ll = thick_max; 
                end
                if (ll>0)
                    if r>c, aar = -aa;
                    else aar = aa;
                    end
                    offs = aar*[-rtheta(r,c) sqrt(1-rtheta(r,c)^2)];
                    if sign(GCJs(Tid,r,c))== +1
                        col = [0.75 0 0.35];
                    elseif sign(GCJs(Tid,r,c))== -1
                        col = 'b';
                    end
                    arrow([crdsf(c,1)+offs(1) crdsf(c,2)+offs(2)],[crdsf(r,1)+offs(1) crdsf(r,2)+offs(2)],'Color',col,'Width',ll,'Length',15,'TipAngle',30);
                end
            end
        end
    end
    axis([0 ymax 0 ymax])
    axis xy
    set(gca,'xticklabel',{[]}, 'yticklabel', {[]});
    title(['GC maps for Tone f='  num2str(2.^((Tid-1)/2)*4) 'kHz @ ' mode ' mode, FDR= ' num2str(alpha) ', N= ' num2str(N)])
    %saveas(gcf,['GCnetworkMapsTone' num2str(2.^((Tid-1)/2)*4) 'mode' mode 'FDR' num2str(alpha) 'N' num2str(N) 'Mod9CV.fig']);
    %close('gcf')
end


