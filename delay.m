function [signalDelayed] = delay(signal,modelOrder,delaySize)
%DELAY Summary of this function goes here
%   Detailed explanation goes here

%% Initialise parameters for use in loop
N = length(signal);
signalDelayed = zeros(modelOrder,N);

%% Loop over model order
for i = 1:modelOrder
    % each row contains the original signal padded by an incremental zeros
    signalDelayed(i,:) = [zeros(1,(i-1)+delaySize) signal(1:N-(i-1)-delaySize)];
end
