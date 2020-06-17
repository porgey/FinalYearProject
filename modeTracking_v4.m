function [modes] = modeTracking_v4(frequency,powerSpectrum,modeStruct,plotFlag)
%modeTracking_v3 Track a specified spectral mode in a spectrogram
%   [mode] = modeTracking_v3(frequency,powerSpectrum,modeStart,globalMin,golbalMax,localBandwidth,time)
%   returns the tracked frequency mode extracted from the power spectrum
%   'powerSpectrum' and its corresponding frequency measurement vector
%   'frequency', using the inout parameters set by the user.
%
%   Inputs:
%       frequency: Frequencies at which power spectrum is computed.
%
%       powerSpectrum: Power spectrum of signal computed by spectrogram.
%
%       modeStruct:
%
%           modeStart: Approximate start of specified frequency mode.
%
%           globalMin: Approximate overall minimum of specified mode
%
%           globalMax: Approximate overall maximum of specified mode
%
%           localBandwidth: Two-sided bandwidth around previous instantaneous
%           frequency in which the current instantaneous frequency will be
%           searched for.
%
%           time - time series of spectrogram, used for plotting
%           
%           solStr - String used for plotting correct Sol
%
%           axisStr - String used for plotting correct axis
%
%       plotFlag - Flag used to specify whether to plot locally-thresholded modes
%
%   Outputs:
%       mode: Time series of the tracked frequency mode.

% get variables applicable to all modes from modeStruct
time = modeStruct.time;
solStr = modeStruct.solStr;
axisStr = modeStruct.axisStr;

% number of modes
nModes = length(modeStruct.modeStarts);
% number of time samples
[~,nTime] = size(powerSpectrum);

%initialise output modes array
modes = zeros(nTime,nModes); 

if plotFlag
    figure;
    sgtitle([solStr ' ' axisStr ' - Locally Thresholded Modes'],'FontWeight','bold','FontSize',14);
end

%% loop over every specified mode
for i = 1:nModes
    % get mode-specific variables from modeStruct
    modeStart = modeStruct.modeStarts(i);
    globalMin = modeStruct.globalMins(i);
    globalMax = modeStruct.globalMaxs(i);
    localBandwidth = modeStruct.localBandwidths(i);

    % define global frequency search area
    fSearchGlobal = (frequency>=globalMin & frequency<=globalMax);
    
    % % Localised thresholding
    % using the 75% quantile
    quantileVal = 0.75;  
    % initialise thresholded spectrum
    powerSpectrum_thresholded = powerSpectrum;
    % find smallest value in spectrogram
    minVal = min(min(powerSpectrum));
    % set power values outside flobal search area to minVal
    powerSpectrum_thresholded(~fSearchGlobal,:) = minVal;
    
    for j = 1:nTime
        % at each time instant, find the upper quantile
        upperQuantile = quantile(powerSpectrum(fSearchGlobal,j),quantileVal);
        
        % set values below the upper quantile to minVal
        powerSpectrum_thresholded(powerSpectrum(:,j) < upperQuantile,j) = minVal;
    end
    
    % % Iterative mode tracking
    % initialise local search area as mode centre and preallocate mode
    fSearchLocalCentre = modeStart;    
    mode = zeros(1,nTime);
    
    for j = 1:nTime
        % obtain local search area around previous mode instantaneous frequency
        fSearchLocal = (frequency>=(fSearchLocalCentre-localBandwidth/2) & frequency<=(fSearchLocalCentre+localBandwidth/2));
        
        % clip local search area by global search area
        fSearchLocal_clipped = fSearchLocal & fSearchGlobal;
        
        % compute median frequency in clipped local search area
        fSearchLocalCentre = medfreq(powerSpectrum(fSearchLocal_clipped,j),frequency(fSearchLocal_clipped));
        
        % set mode at this time instance as the median frequency
        mode(j) = fSearchLocalCentre;
    end
    
    % set current column of output modes array to the current extracted mode 
    modes(:,i) = mode;

    % Plot thresholded spectrum and extracted mode if specified
    if plotFlag
%         figure;
        subplot(1,2,i);
        tStart = time(1);
        tEnd = time(end);
        colormap jet;
        hold on;
        box on;
        surf(time,frequency,log10(sqrt(powerSpectrum_thresholded)),'edgecolor','none');
        c = colorbar;
        view(0,90);
        caxis([-9 -7.75]);
        xlim([tStart tEnd]);
        ylim([globalMin-(1/8)*(globalMax-globalMin) globalMax+(1/8)*(globalMax-globalMin)]);
        shading(gca,'flat');
        title([num2str(modeStart,'%.1f') 'Hz Mode'],'FontSize',14);
%         title({[solStr ' ' axisStr ],['Locally Thresholded ' num2str(modeStart,'%.1f') 'Hz Mode']},'FontSize',14);
        xlabel('Time (LMST)','FontSize',12);
        ylabel('Frequency (Hz)','FontSize',12);
        ylabel(c,'ASD (log(m/s^2/\surd{Hz}))','FontSize',12);
        datetick('x','HH:MM:SS', 'keeplimits');
        plot(time,mode,'Color','black');
        hold off;
    end
    
end

% Plot extracted mode separately
% if plotFlag
%     figure;
%     for i = 1:nModes
%         subplot(1,nModes,i);
%         plot(time,modes(:,i));
%         ylim([modeStruct.globalMins-0.5*(modeStruct.globalMaxs-modeStruct.globalMins)...
%             modeStruct.globalMaxs+0.5*(modeStruct.globalMaxs-modeStruct.globalMins)]);
%         xlim([tStart tEnd]);
% 
%         title(['Extracted ' num2str(modeStruct.modeStarts(i),'%.1f') ' Hz Mode'],'FontSize',14);
%         xlabel('Time (LMST)','FontSize',12);
%         ylabel('Frequency (Hz)','FontSize',12);
%         datetick('x','HH:MM:SS', 'keeplimits');
%     end
% end

end
