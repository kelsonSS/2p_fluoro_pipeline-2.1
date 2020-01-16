function [IGClow,IGChigh] = FDRcontrolBHv4(Devtotal,DkSmth,alpha,Md,wftotal,Whc)
%%% BH rule : Benjamini-Hochberg rule 
% Mod v3: Apr 6th 17: Adapt to static case with input 
% Devtotal of size [Ncells x Ncells x Nmodes/Samples] & 
% wftotal cell of size : [Mf x 1]
% No Mu intercept parameter

[Ncells,~,N] = size(Devtotal);
Mhc = numel(Whc);
[Mf,~] = size(wftotal{1}); Mhcs = Mf - (Ncells-1)*Mhc;

Nl = Ncells*(Ncells-1); % Number of links/OR multiple tests
pth = alpha*(1:Nl)/(Nl*log(Nl)); 
alpham = alpha*(Nl+1)/(2*Nl*log(Nl)); % mean FDR by BH rule

IGC = zeros(size(Devtotal)); IGClow = IGC; IGChigh = IGC;
% Compute p-values
Pv = chi2cdf(Devtotal,Md,'upper');

for t = 1:N
    Pvt = Pv(:,:,t);
    [Pvtsort, Indsrt] = sort(Pvt(:));
    Pvtsortx = Pvtsort(1:Nl);
    cnt = 0;
    while Pvtsortx(cnt+1) < pth(cnt+1)
        cnt = cnt+1;
        if cnt == Nl
            break
        end
    end
    siglind = Indsrt(1:cnt);
    cols = ceil(siglind/Ncells);
    rows = siglind - (cols-1)*Ncells;
    for ee = 1:numel(siglind)
        IGC(rows(ee),cols(ee),t) = 1 - alpham - ncx2cdf( chi2inv(1-alpham,Md),Md, DkSmth(rows(ee),cols(ee),t) );
        whati = wftotal{rows(ee),t};
        CellC = 1:Ncells; CellC(rows(ee))=[];
        jji = find(CellC == cols(ee));

        wseff = whati((jji-1)*Mhc+Mhcs+1:jji*Mhc+Mhcs).*Whc';
        if sign(sum(wseff)) ~=0
            %IGCs(rows(ee),cols(ee),t) = sign(sum(wseff))*IGC(rows(ee),cols(ee),t);
            % Modified on June 1st 17: For more effective and better
            % detection of nature of GC links, take effective nature of
            % links in lower and higher latencies separately. (the lower
            % latencies are better indicators of nature of GC links.)
            IGClow(rows(ee),cols(ee),t) = sign(sum(wseff(1:2)))*IGC(rows(ee),cols(ee),t);
            IGChigh(rows(ee),cols(ee),t) = sign(sum(wseff(3:end)))*IGC(rows(ee),cols(ee),t);
        else
            IGClow(rows(ee),cols(ee),t) = IGC(rows(ee),cols(ee),t); % No sgn change if over regularized
            IGChigh(rows(ee),cols(ee),t) = IGC(rows(ee),cols(ee),t); % No sgn change if over regularized
        end
    end
end
