% Script to measure ambient spectral power distribution using the PR-670
% scanning whatever is in its field of view.
% No display colors are presented — screen stays black/closed.
%
% 10/25/2016 rte  Original version with gun measurements
% 2026       Modified: ambient/scene scanning only, no PTB color display

close all;
clearvars;
commandwindow;

%% --- User Input ---
infoinput = inputdlg({'Measurement Number ='}, ':)', [1 50], {'0'});
nMeasure  = cell2mat(infoinput);

nScansInput = inputdlg({'Number of scans to take ='}, ':)', [1 50], {'1'});
nScans = str2double(nScansInput{1});

delayBetweenScans = 2;   % seconds between each scan

MonitorName = 'OmniStudioX31-5_Gen2_White_Device01'
% MonitorName = 'OmniStudioX31-2_Gen2_c2x3_White_02'
% MonitorName = 'OmniStudioX31-2_Gen2_c2x3_White_01'
% MonitorName = 'OmniStudioX31-2_Gen1_c2x3_White_01'
% MonitorName = 'SurfacePro8_c2x3_White_01'
% MonitorName = 'SurfacePro8_Cal_White';

%% --- Wavelength setup (Brainard convention) ---
nWaves      = 81;
wavelengths = linspace(380, 380 + 5*(nWaves-1), nWaves);

%% --- Initialize PR-670 ---
CMCheckInit(5);   % meterType 5 = PR-670

%% --- Prepare output matrix ---
% Columns: wavelength | spd_1 | qual_1 | spd_2 | qual_2 | ...
nCols      = 1 + 2*nScans;
outputData = zeros(nWaves, nCols);
outputData(:,1) = wavelengths';

fprintf('\nStarting ambient scans. Point PR-670 at target.\n');
fprintf('Taking %d scan(s) with %d s delay between each.\n\n', nScans, delayBetweenScans);

%% --- Scan loop ---
for s = 1:nScans
    fprintf('Scan %d of %d ... ', s, nScans);
    Beeper;

    [spd, qual] = PR670measspd;          % measure ambient SPD
    qual_col    = repmat(qual, nWaves, 1);

    outputData(:, 2*s)     = spd;
    outputData(:, 2*s + 1) = qual_col;

    fprintf('done. Quality code: %d\n', qual);

    if s < nScans
        WaitSecs(delayBetweenScans);
    end
end

fprintf('\nAll scans complete.\n');

%% --- Save to text file ---
OutFileName = strcat(MonitorName, '_Spectrum', nMeasure, '.txt');
ResultFileName = strcat(MonitorName, '_Result', nMeasure, '.csv')
fileId      = fopen(OutFileName, 'w');

tabsOut = repmat('%f\t', 1, nCols-1);
for row = 1:nWaves
    fprintf(fileId, [tabsOut, '%f\n'], outputData(row,:));
end
fclose(fileId);
fprintf('Data saved to: %s\n', OutFileName);

%% --- Plot all scans ---
figure; hold on;
colors = lines(nScans);
xx = outputData(:,1);
for s = 1:nScans
    plot(xx, outputData(:, 2*s), 'Color', colors(s,:), ...
         'DisplayName', sprintf('Scan %d (qual=%d)', s, outputData(1, 2*s+1)));
end
hold off;
xlabel('Wavelength (nm)');
ylabel('Spectral Power');
title('PR-670 Ambient SPD Measurements');
legend show;
grid on;
