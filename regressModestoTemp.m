function [temp_est,polyCoeff,correlCoeffs] = regressModestoTemp(modes,modeTime,windTempData,tempTime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% initialise parameters
regressionOrder = 1;
[~,nModes] = size(modes);
[~,nWindTemps] = size(windTempData);

% preallocate variables for use in loop
modeInterp = cell(1,nWindTemps);
polyCoeff = cell(1,nWindTemps);
temp_est = cell(1,nWindTemps);
correlCoeffs = cell(1,nWindTemps);

% loop over number of wind temperatures
for i = 1:nWindTemps
    % loop over number of modes
    for j = 1:nModes
        % use linear interpolation to align mode vector to temperature
        modeInterp{i}(:,j) = interp1(modeTime,modes(:,j),tempTime);
        
        % generate best fit polynomial coefficients from regression
        polyCoeff{i}(:,j) = polyfit(modeInterp{i}(:,j),windTempData(:,i),regressionOrder);
        
        % synthesise temperature using polynomial refression coefficients 
        temp_est{i}(:,j) = polyval(polyCoeff{i}(:,j),modeInterp{i}(:,j));
        
        % compute corelation coefficients matrix and take off-axis element
        [correlCoeffsMtx,~] = corrcoef(modeInterp{i}(:,j),windTempData(:,i));
        correlCoeffs{i}(j) = correlCoeffsMtx(1,2);
    end
end

% % Linear regression of mode onto temperature
% figure;
% hold on;
% box on;
% scatter(modeInterp{1},windTempData,1);
% title('Mode-Temperature Linear Regression','FontSize',14);
% xlabel('Mode Instantaneous Frequency (Hz)','FontSize',12);
% ylabel(['TWINS Temperature (' char(176) 'C)'],'FontSize',12);
% l = lsline();
% l.Color = 'red';
% l.LineWidth = 2;
% l.LineStyle = '--';
% legend('Scatter Data','Linear Fit','FontSize',11);
% xlim([-inf 7]);
% % plot(
% % plot(modeInterp{1},[windTempData ones(length(windTempData),1)]*p{i}');
% hold off;

end

