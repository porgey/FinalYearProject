function [cleanedSignal] = removeTicker_v3(originalSignal,tickerFreqs,fs,stepsize,filterOrder)

%% Initialise algorithm parameters
N = length(originalSignal);
t = 0:1/fs:((N-1)*(1/fs));
amplitude = 1;
nTickerFreqs = length(tickerFreqs);

%% Loop over number of specified ticker frequencies
for i = 1:nTickerFreqs
    % print to command line
    disp(['Removing ticker frequncy at ' num2str(tickerFreqs(i)) 'Hz']);
    
    % create secondary reference sinusoid epsilon and delay vector u
    secondaryNoise = amplitude*sin(2*pi*tickerFreqs(i)*t);
    
    % create  delay vector u using secondary noise
    secondaryNoiseDelayed = delay(secondaryNoise,filterOrder,1);
    
    % estimate primary noise eta using LMS auxiliary function
    [primaryNoiseEst,~,~] = LMS(originalSignal,secondaryNoiseDelayed,stepsize);
    
    % bandpass estimated primary noise to remove additional components
    primaryNoiseEst_filt = bandpass(primaryNoiseEst,[tickerFreqs(i)-0.25 tickerFreqs(i)+0.25],fs,'StopbandAttenuation',80);
    
    % subtract estimated primary noise from original signal
    originalSignal = originalSignal - primaryNoiseEst_filt';
end

% return signal with all specified ticker frequencies removed
cleanedSignal = originalSignal;

end