function [primaryNoiseEst,error,weights] = LMS(originalSignal,secondaryNoiseDelayed,stepSize)

% get filterOrder and signal length from size of delayed signal
[filterOrder,N] = size(secondaryNoiseDelayed);

% preallocate variables for use in loop
weights = zeros(filterOrder,N+1);
primaryNoiseEst = zeros(1,N);
error = zeros(1,N);

% loop over number of time samples in signal
for i = 1:N
    % estimate primary noise as weighted sum of secondary noise samples
    primaryNoiseEst(i) = weights(:,i)' * secondaryNoiseDelayed(:,i);
    
    % compute error signal
    error(i) = originalSignal(i) - primaryNoiseEst(i);
    
    % update filter weights
    weights(:,i+1) = weights(:,i) + stepSize*error(i)*secondaryNoiseDelayed(:,i);
end