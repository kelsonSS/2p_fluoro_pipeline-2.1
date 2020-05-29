function [gammaOpt,Jcost] = CrossValidModTrial2(Resp,gammacv,Whcv,Whscv,LLcv,LRcv)
%%% Leave-One-Out Cross-Validation: Multi-trial mode

%%% Cross-validation Params
% gammacv = [0:0.1:3];
% Whcv = {Whc};
% Whscv = {Whcs};
% LLcv = 100 %LL
% LRcv = LR

[N,R,Ncells] = size(Resp);
Jcost = zeros(Ncells,numel(gammacv),numel(Whcv));

for ii = 1:numel(Whcv)
    Whi = Whcv{ii}; Mhc = numel(Whi); Lhc = sum(Whi); 
    % Prepare Self-Hist Conditions
    Whsi = Whscv{ii}; Mhcs = numel(Whsi); Lhcs = sum(Whsi);
    Lm = max([Lhc,Lhcs]);
    Np = N - Lm;
    %%% Form Cross-History Covariate Matrices
    Xcz = zeros(Np,Mhc,Ncells,R);
    for r = 1:R
        Xcz(:,:,:,r) = FormHistMatrixv2(squeeze(Resp(:,r,:)),Whi,Lm); 
    end
    
    %%% Form Observation arrays for all cells
    Robs = Resp(Lm+1:end,:,:);
    % Zero-mean the observation vectors
    for cc = 1:Ncells
        Robs(:,:,cc) = Robs(:,:,cc) - ones(Np,1)*mean(Robs(:,:,cc));
    end

    for ct = 1:Ncells
        
        N_indto = num2str(ct);
        x_time = toc;
        tic
        line_length = fprintf(['Cross-validation on cell ' N_indto '...' num2str(floor(x_time)) ' seconds elapsed ']);        
        cellC = 1:Ncells; cellC(ct) = [];        
        
        %%% Form SELF-History Covriates
        Xcsz = zeros(Np,Mhcs,R);
        for r = 1:R
            Xcsz(:,:,r) = FormHistMatrixv2(squeeze(Resp(:,r,ct)),Whsi,Lm); 
        end
        
        %%% Form Effective response vector
        robsz = Robs(:,:,ct);
        
        % Form Full model design/covariate matrix Xf
        Xfc = Xcsz;
        for cc = 1:Ncells-1
            Xfc = cat(2,Xfc,squeeze(Xcz(:,:,cellC(cc),:)));
        end
        Xf = Xfc; 
        
        %%% Prepare for Cross-Validation
        sig2test = zeros(R,numel(gammacv));
        %%% Training stage on all but one repetition (r)
        for r = 1:R
            rC = 1:R ; rC(r) = [];
           % disp(['CV @ repetiton # ' num2str(r) '...'])

            % Form Testing Data
            Xtest = Xf(:,:,r);
            rtest = robsz(:,r); rntest = rtest/sqrt(var(rtest));
            %%% Form Training Data
            % Form Effective standardized Full Model Covariate Matrix
            [Xntrain,Dntrain] = FormFullDesignv3(Xcz(:,:,:,rC),Xcsz(:,:,rC),ct);
            rtrain = robsz(:,rC); rtrain = rtrain(:);
            rntrain = rtrain/sqrt(var(rtrain));
            
            for jj = 1:numel(gammacv)
                gammaj = gammacv(jj);
                %disp(['CV @ gamma = ' num2str(gammaj) '...'])
                % Training Stage: Filtering
                [wntrain,~] = StaticEstimCV(Xntrain,rntrain,gammaj,LLcv,LRcv);            
                wtrain = Dntrain\wntrain;
                %%% Test stage: Compute the variance for the given repetition
                sig2test(r,jj) = (norm(rntest - Xtest*wtrain)^2)/Np;
            end
            
        end
        Jcost(ct,:,ii) = nanmean(sig2test);
        fprintf(repmat('\b',1,line_length))
    end
end

gammaOpt = zeros(Ncells,1);
for ct = 1:Ncells
    LLct = Jcost(ct,:,1);
    [indgamma] = find(LLct == min(min(LLct)));
    if numel(indgamma) > 1
        indgamma = max(indgamma);
    end
    gammaOpt(ct) = gammacv(indgamma);
end

end

