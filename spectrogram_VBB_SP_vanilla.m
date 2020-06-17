close all;
clear;

addpath(genpath('/Users/George/OneDrive - Imperial College London/Fourth Year/Final Year Project/Code'));
addpath(genpath('/Users/George/Google Drive/code'));
load('/Users/George/OneDrive - Imperial College London/Fourth Year/Final Year Project/Code/FYPdata/Sol 99/eStruct.mat');

%% Download Data
% % load Sol194to195_VBB_ZNE.mat 
% Sol1 = 99; %239
% [datetime_sol_begins, ~] = SPsol2utcstartend(Sol1);
% i = 1;
% 
% %Time of data
% e(i).t_event = datetime_sol_begins;%datetime(2019,10,26,06,58,58);
% t1 = 60*60*1; %time before event
% t2 = 60*60*24; %time after event
% 
% %Load and correct 
% e = SPeventdata_dev2(e, t1, t2, 100);
% 
% method = 1;
% [e] = CleanMultipleSeisForUse(e,method);
% e = SPraw2velaccel_dev(e);
% 
% %Load Met data
% e = metdata_dev(e, t1, t2);
% %Align Data and convert to same axis
% e = Align_DataVectors(e);
% e = Met2SEISTime(e);
% 
% sound(sin(1:3000));

%% Extract Data

sensor = 'VBB';

axisExtract = 3;

VBB_data = e.vbbzne(axisExtract).a;     % 1:Z, 2:N, 3:E
SP_data = e.spzne(axisExtract).a;
t_VBB = e.vbb(1).t;
t_SP = e.sp(axisExtract).t;

try
    if isequal(t_VBB,t_SP)
        t_SP = e.t_VBB_LMST;
    else
        t_SP = e.t_SP_LMST;
    end
    t_VBB = e.t_VBB_LMST;
    sol = e.sol;
catch
    [sol,t_VBB_LMST] = SPutc2lmst(t_VBB);
    t_VBB = datetime(datenum(t_VBB_LMST),'ConvertFrom','datenum');
    [~,t_SP_LMST] = SPutc2lmst(t_SP);
    t_SP = datetime(datenum(t_SP_LMST),'ConvertFrom','datenum');
    
    e.t_VBB_LMST = t_VBB;
    e.t_SP_LMST = t_SP;
    e.sol = sol;
end

[solStr,axisStr] = generateSolStr_v2(sol,axisExtract,sensor);

fsVBB = e.vbbsr_max;
fsSP = e.spsr_max;

switch sensor
    case 'VBB'
        SEIS_data = VBB_data;
        t_SEIS = t_VBB;
    case 'SP'
        SEIS_data = SP_data;
        t_SEIS = t_SP; 
end

N = length(SEIS_data);


%% Spectrogram
tinterval = 100;
NaverageSpec = 1.5;
div = 1;
div2 = 1;
sampleNumber = tinterval*fsVBB;
w  = hann(floor(sampleNumber/div)); % hanning window

[~,f,t,pZ] = (spectrogram(SEIS_data, w, ...
    floor(sampleNumber/NaverageSpec/div)/div2, ...
    floor(sampleNumber/NaverageSpec)/div2, ...
    fsVBB,'yaxis'));

pZ_logRoot = log10(sqrt(pZ));

time = t_SEIS(1)+seconds(t);

tStart = time(1);
tEnd = time(end);

%% Plotting
figure;
% subplot(2,1,1);
colormap jet;
box on;
hold on;
surf(time,f,pZ_logRoot,'edgecolor','none');
c = colorbar;
view(0,90);
caxis([-9 -6.5]);
xlim([tStart tEnd]);
if strcmp(sensor,'VBB')
    ylim([0 10]);
elseif strcmp(sensor,'SP')
    ylim([10 50]);

end
shading(gca,'flat');

title([solStr ' ' axisStr ' ' char(8211) ' Axial Spectrogram'],'FontSize',14);
xlabel('Time (LMST)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
ylabel(c,'ASD (log(m/s^2/\surd{Hz}))','FontSize',12);
datetick('x','HH:MM:SS', 'keeplimits');
hold off;


% %% Extract VBB Data
% 
% sensor = 'VBB';
% 
% axisExtract = 3;
% 
% VBB_data = e.vbbzne(axisExtract).a;     % 1:Z, 2:N, 3:E
% SP_data = e.spzne(axisExtract).a;
% t_VBB = e.vbb(1).t;
% t_SP = e.sp(axisExtract).t;
% 
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
% 
% [solStr,axisStr] = generateSolStr_v2(sol,axisExtract,sensor);
% 
% fsVBB = e.vbbsr_max;
% fsSP = e.spsr_max;
% 
% switch sensor
%     case 'VBB'
%         SEIS_data = VBB_data;
%         t_SEIS = t_VBB;
%     case 'SP'
%         SEIS_data = SP_data;
%         t_SEIS = t_SP; 
% end
% 
% N = length(SEIS_data);


% %% Spectrogram
% tinterval = 100;
% NaverageSpec = 1.5;
% div = 1;
% div2 = 1;
% sampleNumber = tinterval*fsVBB;
% w  = hann(floor(sampleNumber/div)); % hanning window
% 
% [~,f,t,pZ] = (spectrogram(SEIS_data, w, ...
%     floor(sampleNumber/NaverageSpec/div)/div2, ...
%     floor(sampleNumber/NaverageSpec)/div2, ...
%     fsVBB,'yaxis'));
% 
% pZ_logRoot = log10(sqrt(pZ));
% 
% time = t_SEIS(1)+seconds(t);
% 
% tStart = time(1);
% tEnd = time(end);

% %% Plotting
% % figure;
% subplot(2,1,2);
% colormap jet;
% box on;
% hold on;
% surf(time,f,pZ_logRoot,'edgecolor','none');
% c = colorbar;
% view(0,90);
% caxis([-9 -6.5]);
% xlim([tStart tEnd]);
% if strcmp(sensor,'VBB')
%     ylim([0 10]);
% elseif strcmp(sensor,'SP')
%     ylim([10 50]);
% 
% end
% shading(gca,'flat');
% 
% title([solStr ' ' axisStr ' ' char(8211) ' Axial Spectrogram'],'FontSize',14);
% xlabel('Time (LMST)','FontSize',12);
% ylabel('Frequency (Hz)','FontSize',12);
% ylabel(c,'ASD (log(m/s^2/\surd{Hz}))','FontSize',12);
% datetick('x','HH:MM:SS', 'keeplimits');
% hold off;

% %% Plot 4Hz Mode
% figure;
% colormap jet;
% box on;
% hold on;
% surf(time,f,pZ_logRoot,'edgecolor','none');
% c = colorbar;
% view(0,90);
% caxis([-9 -6.5]);
% xlim([tStart tEnd]);
% ylim([3.5 4.5]);
% shading(gca,'flat');
% 
% title([solStr ' ' axisStr ' Axial Spectrogram ' char(8211) ' 4.0Hz Mode'],'FontSize',14);
% xlabel('Time (LMST)','FontSize',12);
% ylabel('Frequency (Hz)','FontSize',12);
% ylabel(c,'ASD (log(m/s^2/\surd{Hz}))','FontSize',12);
% datetick('x','HH:MM:SS', 'keeplimits');
% hold off;
