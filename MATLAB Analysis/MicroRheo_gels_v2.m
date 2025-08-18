% Microrheology of gels
% script to process data collected using MicroManager real-time acquisition
% Author: Tania Mendonca
%--------------------------------------------------------------------------

clear all
close all

% Load data
[DataDir] = uigetdir;
[~,name] = fileparts(DataDir);
data = bead_loadData(name);
[Path,name] = fileparts(DataDir);
RFolder = strcat(Path,'\MSD\');    % create output directory

% Compute PSD to check noise
data = bead_preProcessCentres(data);
data = bead_normMSD(data, 'doPlots', false);
data = bead_PSD(data);

%% set variables interactively
% Set Filter Noise to 'n' if want to bypass de-noising
% Noise peaks currently hard-coded to be replaced with noise peaks from
% your data
prompt = {'Bead Radius (m)', 'Sample Temperature (^oC)',...
    'Filter Noise? (y/n):','Frequencies to Filter'};
dlgtitle = 'Inputs';
dims = [1 100];
definput = {'3e-6','37','y'...
    '19.9; 66.3; 332'};        
Inputs = (inputdlg(prompt,dlgtitle,dims,definput));

a = str2num(cell2mat(Inputs(1)));
Temp = str2num(cell2mat(Inputs(2)));
filterCondition = Inputs(3);
bandstopfrequencies = str2num(cell2mat(Inputs(4)));

% Display Raw data plots
figure;                    % Scatter plot
tiledlayout(3,2)
nexttile
plot(data.pro.xCentresM,data.pro.yCentresM); axis equal;
xlabel('X [m]'); ylabel('Y [m]');
title('Raw Scatter')

nexttile                   % XYtime trajectories
plot(data.pro.timeVecMs,data.pro.xCentresM);
hold on; plot(data.pro.timeVecMs,data.pro.yCentresM);
xlabel('Time [ms]'); ylabel('Displacement [m]');
legend({'X','Y'}); title('Raw Trajectories')

nexttile ([2 2])           % MSD from raw data (computed using Will's code)
loglog(data.pro.aMSDnorm(:,1),data.pro.aMSDnorm(:,3),'linewidth',2)
hold on; loglog(data.pro.aMSDnorm(:,2),data.pro.aMSDnorm(:,4),...
    'linewidth',2)
xlabel('Lag Time [s]'); ylabel('MSD [\mum^2]');
legend({'X','Y'}); title('Raw MSD')

% De-noise (Set Filter Noise? to 'y' to run)
if strcmp(filterCondition,'y')
    data.pro.xCentresBS = multfilt(data.pro.xCentresM,bandstopfrequencies,...
        max(data.pro.psd(:,1))*2);
    data.pro.yCentresBS = multfilt(data.pro.yCentresM,bandstopfrequencies,...
        max(data.pro.psd(:,1))*2);
    
    % Drift correction (remove long time drift)
    data.pro.xCentresBS = detrend(data.pro.xCentresBS,1);
    data.pro.yCentresBS = detrend(data.pro.yCentresBS,1);
else
    % Drift correction (remove long time drift)
    data.pro.xCentresBS = detrend(data.pro.xCentresM,1);
    data.pro.yCentresBS = detrend(data.pro.yCentresM,1);
end

% Compute MSD from de-noised data
data.opts.UseField = 'CentresBS';
data = bead_normMSD(data,'forceRun',true,'doPlots', false);

% Extract cleaned trajectories and MSD for saving
txy = [data.pro.timeVecMs' data.pro.xCentresBS' data.pro.yCentresBS'];
msdx = [data.pro.aMSDnorm(2:end,1) data.pro.aMSDnorm(2:end,3)];
msdy = [data.pro.aMSDnorm(2:end,2) data.pro.aMSDnorm(2:end,4)];

% Smoothing MSD with a moving mean filter.
% filter window size can be changed; smaller = less smoothing. Default set
% currently to 30. Max is 227 which is size of the MSD data
msdxS = movmean(msdx, 30);  
msdyS = movmean(msdy, 30);

% Display cleaned data plots
figure;
tiledlayout(3,2)
nexttile
plot(data.pro.xCentresBS,data.pro.yCentresBS); axis equal;
xlabel('X [m]'); ylabel('Y [m]');
title('De-noised Scatter')

nexttile
plot(data.pro.timeVecMs,data.pro.xCentresBS);
hold on; plot(data.pro.timeVecMs,data.pro.yCentresBS);
xlabel('Time [ms]'); ylabel('Displacement [m]');
legend({'X','Y'}); title('De-noised Trajectories')

nexttile ([2 2])
loglog(data.pro.aMSDnorm(:,1),data.pro.aMSDnorm(:,3),'linewidth',2)
hold on; loglog(data.pro.aMSDnorm(:,2),data.pro.aMSDnorm(:,4),...
    'linewidth',2)
% Smoothed MSD overlaid on de-noised MSD
loglog(msdxS(:,1),msdxS(:,2),':k','linewidth',2);...
    loglog(msdyS(:,1),msdyS(:,2),':k','linewidth',2)
xlabel('Lag Time [s]'); ylabel('MSD [\mum^2]');
legend({'X','Y','Smoothed X','Smoothed Y'},'Location', 'southeast'); 
title('De-noised MSD')

%% Microrheology

% Estimate G0 from variance
vx = var(data.pro.xCentresBS)    % variance in x
vy = var(data.pro.yCentresBS)    % variance in y

G0X = 1.38E-23 * (273.15 + Temp) ./ (pi * a * vx);
G0Y = 1.38E-23 * (273.15 + Temp) ./ (pi * a * vy);

%% Save Results
% Currently saving cleaned trajectory text file and smoothed MSD 

if ~exist(RFolder)
    mkdir(RFolder)
end

% convert units to s and m^2
txy(:,1) = txy(:,1)*10^-3; txy(:,2:3) = txy(:,2:3)*10^-6;
%msdx(:,2) = msdx(:,2)*10^-12; msdy(:,2) = msdy(:,2)*10^-12;
msdxS(:,2) = msdxS(:,2)*10^-12; msdyS(:,2) = msdyS(:,2)*10^-12;

save(strcat(name,'.txt'),'txy','-ascii','-double','-tabs');
save(strcat(RFolder,name,'msdxSmoothed.txt'),'msdxS','-ascii','-double','-tabs');
save(strcat(RFolder,name,'msdySmoothed.txt'),'msdyS','-ascii','-double','-tabs');