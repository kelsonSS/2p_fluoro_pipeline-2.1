function [toneVec toneVecN levelVec levelVecN] = correctToneLevelFreq(Frqs, Lvls, toneVec, levelVec)
%% toneBoxAnalysis analyses and plots data from individual ToneBoxes.
%NOT NEEDED FOR DATA COLLECTED AFTER 12/8/18, so this is poorly commented.
% Nikolas Francis, 12/2018

%toneVec
if ~isempty(levelVec)
    Fn = str2num(Frqs);
    [j N] = sort([0; Fn],'descend');
    FL=[];
    for i = 1:length(Frqs)
        FL(i) = length(strtrim(Frqs(i,:)));
    end
    [j l] = sort(FL,'descend');
    Fl = [Frqs(l,:); '   0'];
    toneVectemp=toneVec;
    toneVecN=nan(size(toneVec));
    for i = 1:length(Fl)
        idx = strfind(toneVectemp,strtrim(Fl(i,:)));
        for ii = 1:length(idx)
            if ~str2num(Fl(i,:))
                toneVecN(idx(ii))=0;
            else
                iii = find( Fn == str2num(Fl(i,:)) );
                toneVecN(idx(ii))=iii;
            end
            toneVectemp(idx(ii):idx(ii)+(length(strtrim(Fl(i,:)))-1)) = nan;
        end
    end
    toneVecN(isnan(toneVecN))=[];
    toneVec =num2str(toneVecN);
    toneVec= toneVec(find(~isspace(toneVec)));
else
    if contains(toneVec,'.')
        toneVecN = [];
    else
        for i = 1:length(toneVec)
            toneVecN(i) = str2num(toneVec(i));
        end
    end
end

%levelVec
if ~isempty(levelVec)
    Ln = str2num(Lvls);
    [j N] = sort(Ln,'descend');
    LL=[];
    for i = 1:length(Lvls)
        LL(i) = length(strtrim(Lvls(i,:)));
    end
    [j l] = sort(LL,'descend');
    Ll = [Lvls(l,:); '  X'];
    levelVectemp=levelVec;
    levelVecN=nan(size(levelVec));
    for i = 1:length(Ll)
        idx = strfind(levelVectemp,strtrim(Ll(i,:)));
        for ii = 1:length(idx)
            if strcmpi(Ll(i,:),'  X')
                levelVecN(idx(ii))=0;
            else
                iii = find( Ln == str2num(Ll(i,:)) );
                levelVecN(idx(ii))=iii;
            end
            levelVectemp(idx(ii):idx(ii)+(length(strtrim(Ll(i,:)))-1)) = nan;
        end
    end
    levelVecN(isnan(levelVecN))=[];
    levelVec =num2str(levelVecN);
    levelVec= levelVec(find(~isspace(levelVec)));
else
    levelVec = num2str(zeros(size(toneVec)));
    levelVec= levelVec(find(~isspace(levelVec)));
    levelVecN = zeros(size(levelVec));
end


