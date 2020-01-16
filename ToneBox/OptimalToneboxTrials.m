function [optTrials c] = OptimalToneboxTrials(Data, activitywindow,cfact)
%% OptimalToneboxTrials selects the trials based on optimal responding from an ideal observer model (i.e., ROC)
%Nikolas Francis, 12/2018

%Range of analyzed behavioral response rates
RR = 0:.005:1;

%Find performance statisitcs
H=[];
F=[];
trials=[];
for rr = 1:length(RR)
    
    %Loop over stimuli
    responserate = RR(rr);
    
    %Smoothed first response data
    R =  (movmean(Data.H,activitywindow) +  movmean(Data.F,activitywindow) +  movmean(Data.E,activitywindow))/3;
    
    %Find N active periods above response rate.
    r = find(R>responserate);
    trials{rr} = r;
    
    %Vector of breaks in responsieness
    dR = [1; find(diff(r)>1)];
    
    Ridx=0;
    if length(dR)>1
        for iii = 1:length(dR)-1
            
            Ridx = Ridx+1;
            responses = Data.responseVec(r(dR(iii):dR(iii+1)));
            tones = Data.toneVec(r(dR(iii):dR(iii+1)));
            
            if length(Data.target)==1
                T = strfind(tones,num2str(Data.target));
            else
                sil = strfind(tones,'0');
                T = 1:length(tones)-length(sil);
            end
            
            if ~isfield(Data,'nontarget')
                NT = strfind(tones,'0');
            else
                NT = strfind(tones,num2str(Data.target));
            end
            
            %Hit
            H(Ridx,1,rr) = length(strfind(responses,'H'))./(length(T));
            if H(Ridx,1,rr) == 1
                H(Ridx,1,rr) = 1 - (1/(4*length(T)));
            end
            if H(Ridx,1,rr) == 0
                H(Ridx,1,rr) = (1/(4*length(T)));
            end
            
            %False Alarm
            if ~isempty(strfind(responses,'F'))
                F(Ridx,1,rr) = length(strfind(responses,'F'))./(length(NT));
            else
                F(Ridx,1,rr) = length(strfind(responses,'E'))./(length(T));
            end
            if F(Ridx,1,rr) == 0
                F(Ridx,1,rr) = (1/(4*length(NT)));
            end
            if F(Ridx,1,rr) == 1
                F(Ridx,1,rr) = 1-(1/(4*length(NT)));
            end
            
        end
    end
end

%Select the trials where task engagement is maximal. ie. max(d')
ff=roundTo(squeeze((nanmean(F(:,:,1:end-1)))),3);
hh=roundTo(squeeze((nanmean(H(:,:,1:end-1)))),3);
d=norminv(hh)-norminv(ff);
[v m] = max(d);

%Scale the criterion
if cfact == 1
    optTrials = trials{m};
    c = RR(m);
else
    c = RR(m);
    m = find(RR>c*cfact,1,'first');
    c = RR(m);
    optTrials = trials{m};
end
