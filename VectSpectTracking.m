%% Code to download data and examine 
% clear all
% close all

addpath(genpath('/Users/George/OneDrive - Imperial College London/Fourth Year/Final Year Project/Code'));
addpath(genpath('/Users/George/Google Drive/code'));
load('/Users/George/OneDrive - Imperial College London/Fourth Year/Final Year Project/Code/FYPdata/Sol 99/eStruct_Spher_v2.mat');

%% Download and pre-process data
% Sol1 = 71;
% [datetime_sol_begins, ~] = SPsol2utcstartend(Sol1);
% i=1;
% 
% %Time of data
% e(i).t_event = datetime_sol_begins;%datetime(2019,10,26,06,58,58);
% t1 = 60*60*5; %time before event
% t2 = 60*60*24.5*1; %time after event
% 
% %Load and correct 
% e = SPeventdata_dev2(e, t1, t2);
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

%% Calculate vector spectrogram
% 
% i = 1;
% 
% time1 = datetime(2019,02,07,10,00,00);
% 
% 
% % idx1 = floor(length(e(i).sp(1).t));
% idx1 = 1;
% idx2 = floor(length(e(i).sp(1).t));
% 
% % idx1 = find(e(i).sp(1).t>=time1);
% 
% tickerFreqs = 1:1:10;
% 
% SP_Z = e(i).spzne(1).a(idx1(1):idx2(1));
% SP_N = e(i).spzne(2).a(idx1(1):idx2(1));
% SP_E = e(i).spzne(3).a(idx1(1):idx2(1));
% 
% SP_Z = removeTicker_v2(SP_Z,tickerFreqs,e(i).spsrmax);
% SP_N = removeTicker_v2(SP_N,tickerFreqs,e(i).spsrmax);
% SP_E = removeTicker_v2(SP_E,tickerFreqs,e(i).spsrmax);
% 
% tIntervalSpec = 200; % the length of each spectrogram bin in seconds
% angleStepSpec = 15; % the angle search resolution 
% [e(i).spspher.azigramConditioned, e(i).spspher.elegramConditioned, e(i).spspher.powerFinal, e(i).spspher.powerFinalMin, e(i).spspher.t, e(i).spspher.f] = SPazielegram360_gk_copy(e(i).sp(1).t(idx1(1):idx2(1)), SP_E, SP_N, SP_Z, e(i).spsrmax, tIntervalSpec, angleStepSpec);

%% Plot full Azigram

i = 1;
tIntervalSpec = 200; % the length of each spectrogram bin in seconds
angleStepSpec = 15; % the angle search resolution 


timeSP = e(i).sp(1).t;
% [sol,timeSP] = SPutc2lmst(timeSP_UTC);
% timeSP = e.timeSP_LMST;
% sol = e.sol;

try
    timeSP = e.t_VBB_LMST;
    sol = e.sol;
catch
    [sol,timeSP] = SPutc2lmst(timeSP);
    timeSP = datetime(datenum(timeSP),'ConvertFrom','datenum');
    e.t_VBB_LMST = timeSP;
    e.sol = sol;
end

[solStr,~] = generateSolStr(sol);

srSP = e(i).spsrmax;
angleStep = angleStepSpec;
azigramConditioned = e(i).spspher.azigramConditioned;
elegramConditioned = e(i).spspher.elegramConditioned;
powerFinal = e(i).spspher.powerFinal;
t = e(i).spspher.t;
f = e(i).spspher.f;
loglin = 'lin';


figure;
subplot(1,4,[1 2 3]);
box on;
surf(timeSP(1) + seconds(t), f, azigramConditioned);
shading(gca,'flat');
set(gca,'YScale',loglin);
view([0 0 1]);
set(gca,'layer','top');
% colorbar;
% h = colorbar;
title(['Vector Spectrogram (Azimuth) ' solStr],'FontSize',14);
% ylabel(h, 'Azimuth (Degrees)','FontSize', 12, 'FontWeight', 'bold')
xlabel('Time (LMST)','FontSize',12);
ylabel('Frequency (Hz)','FontSize', 12);
datetick('x','HH:MM:SS', 'keeplimits');
colormap hsv;
caxis([0 360]);
ylim([0 srSP/2]);
% ylim([9.5 11.25]);
xlim([timeSP(1) + seconds(t(1)) timeSP(1) + seconds(t(end))]);

subplot(1,4,4);
color_wheel;
title('Azimuth (Degrees)','FontSize',12,'FontWeight','normal');

p.Position = ax.Position;
p.OuterPosition = ax.OuterPosition;

%% Plot zoomed Azigram and mode angular variation

fMin = 3.85;
fMax = 4.2;
filterLengthFreq = 75; %50;
filterLengthAzi = 340; %250; %/19 for 4-4.15Hz
[modeFreq,modeAzi] = trackAzigramShifts(f,azigramConditioned,fMin,fMax,filterLengthFreq,filterLengthAzi);


figure;
subplotPrimary = subplot(3,1,1);
box on;
hold on;
surf(timeSP(1) + seconds(t), f, azigramConditioned);
shading(gca,'flat');
set(gca,'YScale',loglin);
view([0 0 1]);
set(gca,'layer','top');
colorbar;
h = colorbar;
title(['Azigram ' solStr  ', ' num2str(fMin) '-' num2str(fMax) 'Hz Mode'],'FontSize',14);
xlabel('Time (LMST)');
ylabel('Frequency (Hz)','FontSize',12);
ylabel(h, 'Azimuth (Degrees)','FontSize',12)
set(h,'Ticks',0:60:360);
datetick('x','HH:MM:SS', 'keeplimits');
colormap hsv;
caxis([0 360]);
yline(fMin,'--','LineWidth',4,'Color','black');
yline(fMax,'--','LineWidth',4,'Color','black');
xlim([timeSP(1) + seconds(t(1)) timeSP(1) + seconds(t(end))]);
ylim([fMin-((2/3)*(fMax-fMin)) fMax+((2/3)*(fMax-fMin))]);
hold off;


subplotSecondary{1} = subplot(3,1,2);
box on;
hold on;
plot(timeSP(1) + seconds(t),modeFreq,'LineWidth',2);
title('Mode Frequency Variation','FontSize',14);
xlabel('Time (LMST)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
datetick('x','HH:MM:SS', 'keeplimits');
xlim([timeSP(1) + seconds(t(1)) timeSP(1) + seconds(t(end))]);
ylim([min(modeFreq)-0.5*(max(modeFreq)-min(modeFreq)) max(modeFreq)+0.5*(max(modeFreq)-min(modeFreq))]);
% set(subplot2,'yTick',0:30:360);
hold off;


subplotSecondary{2} = subplot(3,1,3);
box on;
hold on;
plot(timeSP(1) + seconds(t),modeAzi,'LineWidth',2);
title('Mode Azimuthal Variation','FontSize',14);
xlabel('Time (LMST)','FontSize',12);
ylabel('Azimuth (Degrees)','FontSize',12);
datetick('x','HH:MM:SS', 'keeplimits');
xlim([timeSP(1) + seconds(t(1)) timeSP(1) + seconds(t(end))]);
ylim([min(modeAzi)-15 max(modeAzi)+15]);
set(subplotSecondary{2},'yTick',0:30:360);
% set(subplotSecondary{2},'yTick',0:60:360);
hold off;

% paste this into the Command Window to align the subplots.
posPrimary = get(subplotPrimary,'Position');
for i = 1:2
    posSecondary = get(subplotSecondary{i},'Position');
    posSecondary(3) = posPrimary(3);
    set(subplotSecondary{i},'Position',posSecondary);
end

