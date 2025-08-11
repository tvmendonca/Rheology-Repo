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
RFolder = strcat(Path,'\MSD2\');

data = bead_preProcessCentres(data);
data = bead_normMSD(data, 'doPlots', false);
data = bead_PSD(data);

%% set variables
prompt = {'Image Pixel Binning', 'Pixel Calibration (um/pix)',...
    'Bead Radius (m)', 'Sample Temperature (^oC)','Filter Noise? (y/n):',...
    'Frequencies to Filter'};
dlgtitle = 'Inputs';
dims = [1 100];
definput = {'1','0.1083','3e-6','37','y'...
    '19.9; 66.3; 332'};
Inputs = (inputdlg(prompt,dlgtitle,dims,definput));

Binning= str2double(cell2mat(Inputs(1))); %% The pixel binning used when acquiring the image. 
xyCalib = str2num(cell2mat(Inputs(2))); %% microns per pixel calibration for each plane
a = str2num(cell2mat(Inputs(3)));
Temp = str2num(cell2mat(Inputs(4)));
filterCondition = Inputs(5);
bandstopfrequencies = str2num(cell2mat(Inputs(6)));

figure;
tiledlayout(3,2)
nexttile
plot(data.pro.xCentresM,data.pro.yCentresM); axis equal;
xlabel('X [\mum]'); ylabel('Y [\mum]');
title('Raw Scatter')

nexttile
plot(data.pro.timeVecMs,data.pro.xCentresM);
hold on; plot(data.pro.timeVecMs,data.pro.yCentresM);
xlabel('Time [ms]'); ylabel('Displacement [\mum]');
legend({'X','Y'}); title('Raw Trajectories')

nexttile ([2 2])
loglog(data.pro.aMSDnorm(:,1),data.pro.aMSDnorm(:,3),'linewidth',2)
hold on; loglog(data.pro.aMSDnorm(:,2),data.pro.aMSDnorm(:,4),...
    'linewidth',2)
xlabel('Lag Time [ms]'); ylabel('MSD');
legend({'X','Y'}); title('Raw MSD')

%% Filter Data (Set 

% De-noise
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

data.opts.UseField = 'CentresBS';
data = bead_normMSD(data,'forceRun',true,'doPlots', false);

txy = [data.pro.timeVecMs' data.pro.xCentresBS' data.pro.yCentresBS'];
msdx = [data.pro.aMSDnorm(2:end,1) data.pro.aMSDnorm(2:end,3)];
msdy = [data.pro.aMSDnorm(2:end,2) data.pro.aMSDnorm(2:end,4)];

msdxS = movmean(msdx, 30);
msdyS = movmean(msdy, 30);

figure;
tiledlayout(3,2)
nexttile
plot(data.pro.xCentresBS,data.pro.yCentresBS); axis equal;
xlabel('X [\mum]'); ylabel('Y [\mum]');
title('De-noised Scatter')

nexttile
plot(data.pro.timeVecMs,data.pro.xCentresBS);
hold on; plot(data.pro.timeVecMs,data.pro.yCentresBS);
xlabel('Time [ms]'); ylabel('Displacement [\mum]');
legend({'X','Y'}); title('De-noised Trajectories')

nexttile ([2 2])
loglog(data.pro.aMSDnorm(:,1),data.pro.aMSDnorm(:,3),'linewidth',2)
hold on; loglog(data.pro.aMSDnorm(:,2),data.pro.aMSDnorm(:,4),...
    'linewidth',2)
loglog(msdxS(:,1),msdxS(:,2),':k','linewidth',2); loglog(msdyS(:,1),...
    msdyS(:,2),':k','linewidth',2)
xlabel('Lag Time [ms]'); ylabel('MSD');
legend({'X','Y','Smoothed X','Smoothed Y'},'Location', 'southeast'); 
title('De-noised MSD')

%% Microrheology

% Estimate G0 from variance
vx = var(data.pro.xCentresBS)
vy = var(data.pro.yCentresBS)

G0X = 1.38E-23 * (273.15 + Temp) ./ (pi * a * vx);
G0Y = 1.38E-23 * (273.15 + Temp) ./ (pi * a * vy);

%% Save Results

if ~exist(RFolder)
    mkdir(RFolder)
end
   
save(strcat(name,'.txt'),'txy','-ascii','-double','-tabs');
save(strcat(RFolder,name,'msdxSmoothed.txt'),'msdxS','-ascii','-double','-tabs');
save(strcat(RFolder,name,'msdySmoothed.txt'),'msdyS','-ascii','-double','-tabs');