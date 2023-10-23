% Time_Series_Normalization

%% Description
% This script requires the input of a list file which includes the
% file name, injection, nadir and tailwater time stamps as ROI.
%
% The ROI are then used to normalize the passage duration by resampling the
% pressure and acceleration magnitude on a scale from 0 to 1.
%
% 0.0 = injection
% 0.5 = nadir pressure
% 1.0 = exit to tailwater
% 
% After completion, the script provides three figures:
% Figure 1: Shows the injection, nadir and tailwater locations in the raw pressure time series.
% Figure 2: Plots an overlay of all pressure time series after normalization.
% Figure 3: Illustrates the overlay of all acceleration magnitude time series after normalization.
%%

close all
clear all

% User-defined sample size of the normalized time series data
N = 100; %

% Select whether to normalize the FBS or BDS data
normalizeFBS = 0; % set this variable to 1 for FBS, or to any other value to will normalize BDS
if normalizeFBS == 1
    % Fish Backpack Sensors (FBS)
    load('C:\Users\Admin\Desktop\2023.10.19_Data\MATLAB\FBS_n6_ROI.mat'); % User location of the ROI list files for FBS
    filePath = 'C:\Users\Admin\Desktop\2023.10.19_Data\MATLAB\FBS_Matlab_Files_n6\'; % User location of the Matlab raw FBS data files
else
    % Barotrauma Detection System Sensors (BDS)
    load('C:\Users\Admin\Desktop\2023.10.19_Data\MATLAB\BDS_n10_ROI.mat'); % User location of the ROI list files for BDS
    filePath = 'C:\Users\Admin\Desktop\2023.10.19_Data\MATLAB\BDS_Matlab_Files_n10\'; % User location of the Matlab raw BDS data files
end

contents = list.name; % requires separate table entry "list" with list.name, list.injection and list.recovery in seconds
num_contents = height(list); % number of files to process

for itFile = 1:size(contents,1) % This iterates over all .mat files with raw sensor data
    fileNameTxt = char(contents(itFile));
    fileNameTxt = fileNameTxt(1:end-4); % remove ending
    fileFull = strcat(filePath,fileNameTxt,'.mat');
    load(fileFull);
    disp(['Importing and transforming file: ',fileNameTxt,' ...']);

timeStart = list.injection(itFile);
timeDiff = list.nadir(itFile) - list.injection(itFile); % time duration from injection until nadir, creates a symmetric plot with the nadir at 50%
timeStop = list.nadir(itFile) + timeDiff;

timeStartPost = list.nadir(itFile);
timeStopPost = list.tailwater(itFile);

if normalizeFBS == 1
    %% Convert FBS files to the same input format as the BDS
    ts = retero.ts; % time in seconds
    P1 = retero.p;
    P2 = P1; P3 = P1; % pressure in mbar
    acc_x = retero.ax; acc_y = retero.ay; acc_z = retero.az;
    %%
else 
end

idxStart = find(ts<=timeStart); % find time stamp in s corresponding to list.deployment starting time
idxStart = idxStart(end);
idxStop = find(ts<timeStop); % find time stamp in s for list.nadir time
idxStop = idxStop(end);

idxStopPost = find(ts<timeStopPost); % find time stamp in s for list.tailwater end time
idxStopPost = idxStopPost(end);

PMean = mean([P1,P2,P3],2); % calculate mean of 3 pressure values at each time step

%% Find true nadir
idxNadir = find(ts<list.nadir(itFile));
idxNadir = idxNadir(end);
nadir_range = 3; % range of values to look for a better nadir
nadir_detect = PMean(idxNadir-nadir_range:idxNadir+nadir_range);
[nadir_value,idx_Nadir_Detect] = min(nadir_detect); % find the lowest value within the range;
idxNadirTrue = idx_Nadir_Detect - nadir_range + idxNadir -1;
list.nadirComputed(itFile) = ts(idxNadirTrue);

%% Create pre and post-nadir dataset which is concatenated
idxStop = idxNadirTrue; % pre nadir data ends with the nadir value
idxStartPost = idxNadirTrue + 1; % post nadir data begins with the first measurement after the nadir
p_cropcat = [PMean(idxStart:idxStop);PMean(idxStartPost:idxStopPost)];
ts_prepost = (1:size(p_cropcat))./(size(p_cropcat,1));
sizePre = idxStop - idxStart + 1;
sizePost = idxStopPost - idxStartPost + 1;
ts_pre = (1:sizePre)/(2*sizePre);
ts_post = (1:sizePost)/(2*sizePost)+0.5;
ts_prepost = [ts_pre,ts_post];
ratioPre = sizePre / (sizePre + sizePost);
ratioPost = sizePost / (sizePre + sizePost);
%%

%% Calculate and resample pressure to fixed passage duration with interval of [0,1]
ts_normalized = ((0:1:N-1)/(N-1))';
pressure_resampled_fixed_interval(:,itFile) = resample(p_cropcat-mean(p_cropcat),ts_prepost,N)+mean(p_cropcat);
pressure_resampled_fixed_interval(1,itFile) = pressure_resampled_fixed_interval(2,itFile); % replace the first value with the exact value from the data, since non-zero entries need to be padded when resampling
pressure_resampled_fixed_interval(end,itFile) = pressure_resampled_fixed_interval(end-1,itFile); % replace the last value with the exact value from the data, since non-zero entries need to be padded when resampling
%%

%% Calculate acceleration magnitude of the cropped data
ax_crop = acc_x(idxStart:idxStop);
ay_crop = acc_y(idxStart:idxStop);
az_crop = acc_z(idxStart:idxStop);
am = sqrt(acc_x.^2 + acc_y.^2 + acc_z.^2);
am_crop = sqrt(ax_crop.^2 + ay_crop.^2 + az_crop.^2);
%%

%% Pre- and post-nadir data are now cropped, resampled and concatenated
acc_mag_cropcat = [am(idxStart:idxStop);am(idxStartPost:idxStopPost)];
acc_mag_fixed_interval(:,itFile) = resample(acc_mag_cropcat-mean(acc_mag_cropcat),ts_prepost,N)+mean(acc_mag_cropcat);
acc_mag_fixed_interval(1,itFile) = acc_mag_fixed_interval(2,itFile); % replace the first value with the exact value from the data, since non-zero entries need to be padded when resampling
acc_mag_fixed_interval(end,itFile) = acc_mag_fixed_interval(end-1,itFile); % replace the last value with the exact value from the data, since non-zero entries need to be padded when resampling 
%%

figure(1); % Figure showing all pressure time series, labelling the injection (green), nadir (black) and tailwater (red)
subplot(num_contents,1,itFile);
plot(ts,PMean,'b','LineWidth',1);
hold on
plot(ts(idxStart),PMean(idxStart),'og');
plot(ts(idxNadirTrue),PMean(idxNadirTrue),'ok');
plot(ts(idxStopPost),PMean(idxStopPost),'or');
ylim([400 1500]);
title(fileNameTxt)
legend('Pressure (mbar)','Injection','Nadir','Tailwater');
xlabel('Time (s)');
ylabel('Total Pressure (mbar)');

figure(2); % Figure showing the pressure time series after normalization from 0 to 1, where nadir is located at 0.5
plot(ts_normalized,pressure_resampled_fixed_interval,'b');
xlabel('Normalized Time (s)');
ylabel('Total Pressure (mbar)');

figure(3); % Figure showing the acceleration magnitude time series after normalization from 0 to 1, where nadir is located at 0.5
plot(ts_normalized,acc_mag_fixed_interval,'r');
xlabel('Normalized Time (s)');
ylabel('Acceleration Magnitude (m/s^2)');
end
