% Calculates bandwidth of the max peak on a curve at an input threshold value.
% Excludes peaks that are discontinous at the boundary.
%   Inputs:     curve - curve to be analyzed (like a tuning curve)
%               peak_threshold - threshold to analyze bandwidth
%               freqs - array of frequency values in kHz
%               show_plot (optional) - enter 1 to produce graph of curve
%   Outputs:    bandwidth - the peak bandwidth. NaN if peak is discontinous
%               bf - the index of the peak
%               Q - the Q factor of the peak
% Zac Bowen Jan2016

function [bandwidth,bf,Q] = peak_analyzer(curve,peak_threshold,kHzVals,show_plot)

if nargin < 4
    show_plot = 0;
end

%% Normalize input curve
% minVal = min(curve);
% tmp = curve - minVal;
% normCurve = tmp / max(tmp);
normCurve = curve / max(curve); %changed norm to 0 to max

%%
octvals = log2(kHzVals/kHzVals(1)); %from kHz to octaves
tonesPerOct = find(octvals==1) - 1; % Finds the tones per octave
[bfamp,bf] = max(normCurve); % Finds peak and peak index
cutoff = bfamp*peak_threshold; % Sets cutoff based on max peak
bandwidth = NaN; % Preallocation - is returned as NaN if peak is on the edge and bandwidth cannot be calculated
Q = NaN;

if bf == 1 || bf == length(curve)
    return
end

% interpolates the curve for higher resolution
res = 100; % increase the resolution of curve 100x
int_curve = interp1(1:length(normCurve),normCurve,1:(1/res):length(normCurve));
int_curve = int_curve';
testline(1:res*(length(normCurve)-1)+1,1) = cutoff; % defines cutoff line

difference = int_curve - testline; %finds differences above/below cutoff
intersects = find(diff(sign(difference)))'; %finds the intersects

% Need to incorporate whether this was the BF peak or not
if int_curve(1) > cutoff && int_curve(end) > cutoff
    intersects = [NaN intersects NaN];
    peakWindowBounds = reshape(intersects,[2,length(intersects)/2]);
%     fprintf('First val and Last val were above threshold\n')
elseif int_curve(1) > cutoff
    intersects = [NaN intersects];
    peakWindowBounds = reshape(intersects,[2,length(intersects)/2]);
%     fprintf('First val was above threshold\n')
elseif int_curve(end) > cutoff
    intersects = [intersects NaN];
    peakWindowBounds = reshape(intersects,[2,length(intersects)/2]);
% %     fprintf('Last val was above threshold\n')
else
    peakWindowBounds = reshape(intersects,[2,length(intersects)/2]);
end

nr_sig_peaks = length(intersects)/2;
if nr_sig_peaks == 0
%     fprintf('No significant peaks\n')
    return
end

[~,int_max] = max(int_curve);
for k = 1:nr_sig_peaks
    if any(int_max == (peakWindowBounds(1,k):peakWindowBounds(2,k)))
        bandwidth = peakWindowBounds(2,k) - peakWindowBounds(1,k);
        bandwidth = (bandwidth - 1)/res;
        peak_kHz_val = kHzVals(1)*2^((bf-1)/tonesPerOct);
        low_kHz_val = kHzVals(1)*2^((peakWindowBounds(1,k)/res)/tonesPerOct);
        high_kHz_val = kHzVals(1)*2^((peakWindowBounds(2,k)/res)/tonesPerOct);
        Q = peak_kHz_val / (high_kHz_val - low_kHz_val);
    end
end

if show_plot
    figure; hold on;
    p = plot(int_curve,'b'); p.Parent.XTickLabel = '';
    plot(testline,'r');
    title(['BF = ' num2str(peak_kHz_val) 'kHz and Bandwidth = ' num2str(bandwidth) ' octaves'])
    ylabel('Normalized response (dF/F %)')
    xlabel('Frequency values (kHz)')
    p.Parent.FontSize = 14;
    p.Parent.XTick = p.Parent.XTick(1:end-1);
    p.Parent.XTickLabel = kHzVals;
end

end