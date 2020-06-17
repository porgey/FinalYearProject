function [modeFreqFiltered,modeAziFiltered] = trackAzigramShifts(frequency,azigramConditioned,fMin,fMax,filterLengthFreq,filterLengthAzi)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% define search area Boolean vector
fSearch = (frequency >= fMin) & (frequency <= fMax);
% find index of first element in fSearch
fSearchLowerIdx = find(fSearch,1,'first');
% define search area in azigram
azigramSearch = azigramConditioned(fSearch,:);
% compute statistical mode at each time sample
modeAzi = mode(azigramSearch);

%% Extract time-frequency mode
nTime = size(azigramConditioned,2);
modeFreq = zeros(1,nTime);

for i = 1:nTime
    % find indices at which the statistical modal azimuth is present
    modeIdxs = find(azigramSearch(:,i)==modeAzi(i));
    
    % estimate instantanous frequency index as the middle index
    freqIdx = modeIdxs(ceil(length(modeIdxs)/2));
    
    % set current instantanous frequency to the corresponding frequency
    modeFreq(i) = frequency(freqIdx+fSearchLowerIdx-1);
end

% median filter the extracted time-frequency mode
modeFreqFiltered = medfilt1(modeFreq,filterLengthFreq,'truncate');

%% Median filter the time-azimuth mode
modeAziFiltered = medfilt1(modeAzi,filterLengthAzi,'truncate');

end