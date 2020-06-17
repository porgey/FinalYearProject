% close all;
clear;

addpath(genpath('/Users/George/OneDrive - Imperial College London/Fourth Year/Final Year Project/Code'));
addpath(genpath('/Users/George/Google Drive/code'));
load('/Users/George/OneDrive - Imperial College London/Fourth Year/Final Year Project/Code/FYPdata/Sol 99/eStruct.mat');

%% Download Data
% % load Sol194to195_VBB_ZNE.mat 
% Sol1 = 361; %239
% [datetime_sol_begins, ~] = SPsol2utcstartend(Sol1);
% i = 1;
% 
% % Time of data
% e(i).t_event = datetime_sol_begins;%datetime(2019,10,26,06,58,58);
% t1 = 60*60*1; %time before event
% t2 = 60*60*24; %time after event
% 
% % Load and correct
% e = SPeventdata_GKdev(e, t1, t2, 100);
% 
% method = 1;
% [e] = CleanMultipleSeisForUse(e,method);
% e = SPraw2velaccel_dev(e);
% 
% % Load Met data 
% e = metdata_GKdev(e, t1, t2);
% % Align Data and convert to same axis
% e = Align_DataVectors(e);
% e = Met2SEISTime(e);
% 
% sound(sin(1:3000));

%% Extract Data

sensor = 'SP';

axisExtract = 3;

% try
%     VBB_data = e.vbbzne(axisExtract).a;     % 1:Z, 2:N, 3:E
%     SP_data = e.spzne(axisExtract).a;
%     t_VBB = e.vbb(axisExtract).t;
%     t_SP = e.sp(axisExtract).t;
% catch
%     SP_data = e.spzne(axisExtract).a;
%     t_SP = e.sp(axisExtract).t;
%     VBB_data = SP_data;
%     t_VBB = t_SP;
% end

    

% try
%     if isequal(t_VBB,t_SP)
%         t_SP = e.t_VBB_LMST;
%     else
%         t_SP = e.t_SP_LMST;
%     end
%     t_VBB = e.t_VBB_LMST;
%     sol = e.sol;
% catch
%     [sol,t_VBB_LMST] = SPutc2lmst(t_VBB);
%     t_VBB = datetime(datenum(t_VBB_LMST),'ConvertFrom','datenum');
%     [~,t_SP_LMST] = SPutc2lmst(t_SP);
%     t_SP = datetime(datenum(t_SP_LMST),'ConvertFrom','datenum');
%     
%     e.t_VBB_LMST = t_VBB;
%     e.t_SP_LMST = t_SP;
%     e.sol = sol;
% end

% [sol,t_SP_LMST] = SPutc2lmst(t_SP);
% t_SP = datetime(datenum(t_SP_LMST),'ConvertFrom','datenum');

% t_VBB = e.t_VBB_LMST;  %these lines assume time has been extraced and converted already
% t_SP = e.t_SP_LMST;
% sol = e.sol;


% try
%     fsVBB = e.vbbsr_max;
%     fsSP = e.spsr_max;
% catch
%     fsSP = e.spsr_max;
%     fsVBB = fsSP;
% end


% VBB_data = highpass(VBB_data,0.25,fsVBB); % remove low-frequency artifacts present in seismic signal (f<~0.25Hz)

% switch sensor
%     case 'VBB'
%         SEIS_data = VBB_data;
%         t_SEIS = t_VBB;
%     case 'SP'
%         SEIS_data = SP_data;
%         t_SEIS = t_SP; 
% end

switch sensor
    case 'VBB'
        try
            SEIS_data = e.vbbzne(axisExtract).a;
            t_SEIS = e.t_VBB_LMST;
            sol = e.sol;
            fs_SEIS = e.vbbsr_max;
        catch
            [sol,t_VBB_LMST] = SPutc2lmst(e.vbb(axisExtract).t);
            t_VBB_LMST = datetime(datenum(t_VBB_LMST),'ConvertFrom','datenum');
            e.t_VBB_LMST = t_VBB_LMST;
            e.sol = sol;
            
            SEIS_data = e.vbbzne(axisExtract).a;
            t_SEIS = e.t_VBB_LMST;
            sol = e.sol;
            fs_SEIS = e.vbbsr_max;
            disp('The "e" struct has new terms. Save to keep the changes.');
        end
    case 'SP'
        try
            SEIS_data = e.spzne(axisExtract).a;
            t_SEIS = e.t_SP_LMST;
            sol = e.sol;
            fs_SEIS = e.spsr_max;
        catch
            try
                %if VBB LMST data already saved and can be used, use that
                if isequal(e.vbb(axisExtract).t,e.sp(axisExtract).t)
                    t_SP_LMST =  e.t_VBB_LMST;
                    e.t_SP_LMST = t_SP_LMST;
                else
                    [sol,t_SP_LMST] = SPutc2lmst(e.sp(axisExtract).t);
                    t_SP_LMST = datetime(datenum(t_SP_LMST),'ConvertFrom','datenum');
                    e.t_SP_LMST = t_SP_LMST;
                    e.sol = sol;
                end
            catch
                % catch errors
                [sol,t_SP_LMST] = SPutc2lmst(e.sp(axisExtract).t);
                t_SP_LMST = datetime(datenum(t_SP_LMST),'ConvertFrom','datenum');
                e.t_SP_LMST = t_SP_LMST;
                e.sol = sol;
            end
            
            SEIS_data = e.spzne(axisExtract).a;
            t_SEIS = e.t_SP_LMST;
            sol = e.sol;
            fs_SEIS = e.spsr_max;
            disp('The "e" struct has new terms. Save to keep the changes.');
        end
end

N = length(SEIS_data);
[solStr,axisStr] = generateSolStr_v2(sol,axisExtract,sensor);


%% Remove clock ticker using ANC
tickerFreqs = [1];
mu = 0.0001;    %0.00001
filterOrder = 10;   %10
plotFlag = 0;
SEIS_data_ticksRemoved = removeTicker_v3(SEIS_data,tickerFreqs,fs_SEIS,mu,filterOrder);


%% Spectrogram
tinterval = 100;
NaverageSpec = 1.5;
div = 1;
div2 = 1;
sampleNumber = tinterval*fs_SEIS;
w  = hann(floor(sampleNumber/div)); % hanning window

[~,~,~,pZ_orig] = (spectrogram(SEIS_data, w, ...
    floor(sampleNumber/NaverageSpec/div)/div2, ...
    floor(sampleNumber/NaverageSpec)/div2, ...
    fs_SEIS,'yaxis'));

[~,f,t,pZ] = (spectrogram(SEIS_data_ticksRemoved, w, ...
    floor(sampleNumber/NaverageSpec/div)/div2, ...
    floor(sampleNumber/NaverageSpec)/div2, ...
    fs_SEIS,'yaxis'));


%% Thresholding
% quantileVal = 0.9;  % Using the 90% quantile
% pZ_thresholded = pZ;
% minVal = min(min(pZ)); % find smallest value in spectrogram
% 
% for i = 1:length(t)
%     upperQuantile = quantile(pZ(:,i),quantileVal);  % at each time instant, find the upper quantile
%     pZ_thresholded(pZ(:,i) < upperQuantile,i) = minVal; % set values below the upper quantile to minVal
% end

pZ_logRoot = log10(sqrt(pZ));
% pZ_logRoot_thresholded = log10(sqrt(pZ_thresholded));

time = t_SEIS(1)+seconds(t);

%% Mode tracking (all at Q=0.75 unless specified)

% modeStruct.modeStarts = [4.05,6.7];
% modeStruct.globalMins = [3.75,5.85];
% modeStruct.globalMaxs = [4.15,7.0];
% modeStruct.localBandwidths = [0.05,0.1];

% modeStruct.modeStarts = [3.35,4.05,6.7,8.4];
% modeStruct.globalMins = [3.05,3.75,5.85,6.8];
% modeStuct.globalMaxs = [3.4,4.15,7.0,8.7];
% modeStruct.localBandwidths = [0.055,0.057,0.09,0.089];

% modeStruct.modeStarts = [25.3];
% modeStruct.globalMins = [24.45];
% modeStruct.globalMaxs = [25.7];
% modeStruct.localBandwidths = [0.17];

modeStruct.modeStarts = [25.3,29.6];
modeStruct.globalMins = [24.45,28.5];
modeStruct.globalMaxs = [25.65,30.9];
modeStruct.localBandwidths = [0.15,0.15];

% modeStruct.modeStarts = [9.8];
% modeStruct.globalMins = [7.0];
% modeStruct.globalMaxs = [10];
% modeStruct.localBandwidths = [0.1];

% modeStruct.modeStarts = [6.7];
% modeStruct.globalMins = [5.85];
% modeStruct.globalMaxs = [7.0];
% modeStruct.localBandwidths = [1.5]; 
modeStruct.time = time;
modeStruct.solStr = solStr;
modeStruct.axisStr = axisStr;
plotFlag = 1;

modes = modeTracking_v4(f,pZ,modeStruct,plotFlag);
[~,nModes] = size(modes);


%% Fitting Modes to Temperature
tStart = time(1);
tEnd = time(end);

try
    temp_time = e.temp_time_LMST;
    scit_time = e.scit_time_LMST;
catch
    [~,temp_time] = SPutc2lmst(e.windtempdata.time);
    temp_time = datetime(datenum(temp_time),'ConvertFrom','datenum');
    e.temp_time_LMST = temp_time;
    
    [~,scit_time] = SPutc2lmst(e.scit.t);
    scit_time = datetime(datenum(scit_time),'ConvertFrom','datenum');
    e.scit_time_LMST = scit_time;
end

%temporary line
% temp_time = scit_time;

idx_start = find(temp_time >= time(1),1,'first');
idx_end = find(temp_time <= time(end),1,'last');
temp_time = temp_time(idx_start:idx_end);

nWindTemps = 1;
windTempData = zeros(length(temp_time),nWindTemps);

windTempData_extracted = e.windtempdata.data_temp;
windTempData(:,1) = windTempData_extracted(idx_start:idx_end);
% windTempData_extracted = e.windtempdata.data_temp_boom1;
% windTempData(:,2) = windTempData_extracted(idx_start:idx_end);
% windTempData_extracted = e.windtempdata.data_temp_boom2;
% windTempData(:,3) = windTempData_extracted(idx_start:idx_end);
% windTempData_extracted = e.scit.d;
% windTempData(:,1) = windTempData_extracted(idx_start:idx_end);

% windTempStr = {'TWINS Composite','TWINS Boom 1','TWINS Boom 2'};
windTempStr = {'TWINS Composite'};
% windTempStr = {'SCIT'};

[temp_est,p,correlCoeffs] = regressModestoTemp(modes,time,windTempData,temp_time);

% % Linear regression of mode onto temperature
% figure;
% hold on;
% scatter(windTempData,modeInterp(:,i));
% plot(modeInterp(:,i),[windTempData ones(length(windTempData),1)]*p(:,i)');
% hold off;


% % Exploring different regression orders
% nOrders = 1;
% temp_est = zeros(length(windTempData),nOrders);
% MSE = zeros(1,3);
% 
% for i = 1:nOrders
%     p = polyfit(modeInterp,windTempData,i);
%     temp_est(:,i) = polyval(p,modeInterp);
%     MSE(i) = (1/length(windTempData))*sum((windTempData-temp_est(:,i)).^2);
% end


%% Plotting
figure;
subplot(2,1,1)
colormap jet;
box on;
hold on;
surf(time,f,log10(sqrt(pZ_orig)),'edgecolor','none');
c = colorbar;
view(0,90);
xlim([tStart tEnd]);
if strcmp(sensor,'VBB')
    ylim([0 10]);
    caxis([-9 -6.5]);
elseif strcmp(sensor,'SP')
    ylim([23 33]);
    caxis([-7.25 -5]);
end
% ylim([modeStruct.globalMins-0.5*(modeStruct.globalMaxs-modeStruct.globalMins)...
%     modeStruct.globalMaxs+0.5*(modeStruct.globalMaxs-modeStruct.globalMins)]);
shading(gca,'flat');
% for i = 1:nModes
%     plot(time,modes(:,i),'color','k','LineWidth',2);
% end
% yline(modeStruct.globalMins,'--','globalMin','LineWidth',4,'Color','black');
% yline(modeStruct.globalMaxs,'--','globalMax','LineWidth',4,'Color','black');
title([solStr ' ' axisStr ' ' char(8211) ' Axial Spectrogram'],'FontSize',14);
xlabel('Time (LMST)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
ylabel(c,'ASD (log(m/s^2/\surd{Hz}))','FontSize',12);
datetick('x','HH:MM:SS', 'keeplimits');
hold off;

subplot(2,1,2);
% subplotPrimary = subplot(nModes+1,1,1);
colormap jet;
box on;
hold on;
surf(time,f,pZ_logRoot,'edgecolor','none');
c = colorbar;
view(0,90);
% caxis([-9 -8]);
xlim([tStart tEnd]);
if strcmp(sensor,'VBB')
    ylim([0 10]);
    caxis([-9 -6.5]);
elseif strcmp(sensor,'SP')
    ylim([23 33]);
    caxis([-7.25 -5]);
end
% ylim([modeStruct.globalMins-0.5*(modeStruct.globalMaxs-modeStruct.globalMins)...
%     modeStruct.globalMaxs+0.5*(modeStruct.globalMaxs-modeStruct.globalMins)]);
shading(gca,'flat');
for i = 1:nModes
    plot(time,modes(:,i),'color','k','LineWidth',2);
end
% yline(modeStruct.globalMins,'--','globalMin','LineWidth',4,'Color','black');
% yline(modeStruct.globalMaxs,'--','globalMax','LineWidth',4,'Color','black');
title([solStr ' ' axisStr ' ' char(8211) ' Tracked Modes'],'FontSize',14);
xlabel('Time (LMST)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
ylabel(c,'ASD (log(m/s^2/\surd{Hz}))','FontSize',12);
datetick('x','HH:MM:SS', 'keeplimits');
hold off;

figure;
sgtitle([solStr ' ' axisStr],'FontSize',14,'FontWeight','bold');
for i = 1:nWindTemps
    subplotSecondary{i} = subplot(nWindTemps,1,i); %#ok<SAGROW>
    box on;
    hold on;
    plot(temp_time,windTempData(:,i),'LineWidth',3);
    legendStr = {windTempStr{i}};
    for j = 1:nModes
        plot(temp_time,temp_est{i}(:,j),'LineWidth',1);
        legendStr{end+1} = [num2str(modeStruct.modeStarts(j),'%.1f') 'Hz ' char(8211) ' r = ' num2str(correlCoeffs{i}(j),'%.3f')]; %#ok<SAGROW>
    end
    xlim([temp_time(1) temp_time(end)]);
    title([windTempStr{i} ' Temperature & Regressed Modes'],'FontSize',14);
%     title('TWINS Temperature & Regressed Modes','FontSize',14);
    xlabel('Time (LMST)','FontSize',12);
    if windTempData(1,i) < 0
        ylabel(['Temperature (' char(176) 'C)'],'FontSize',12);
    elseif windTempData(1,i) > 1000
        ylabel('Temperature (Counts)','FontSize',12);
    else
        ylabel('Temperature (K)','FontSize',12);
    end
    legend(legendStr,'FontSize',11,'Location','northwest');
%     legend('TWINS Temperature','Regressed Mode','FontSize',11,'Location','NorthWest');
    datetick('x','HH:MM:SS', 'keeplimits');
    hold off;
end

% paste this into the Command Window to align the subplots.
% posPrimary = get(subplotPrimary,'Position');
% for i = 1:nModes
%     posSecondary = get(subplotSecondary{i},'Position');
%     posSecondary(3) = posPrimary(3);
%     set(subplotSecondary{i},'Position',posSecondary);
% end

% figure;
% hold on;
% plot(e.vbb1temp.time,e.vbb1temp.raw_data)
% plot(e.scit.time,e.scit.raw_data);
% hold off;

