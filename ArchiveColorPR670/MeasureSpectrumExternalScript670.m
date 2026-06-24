% Script to measure RGB primaries using the PR-670
% GUI-driven: user manually sets each gun, clicks to measure.

close all;
clearvars;
commandwindow;

%% --- User Input ---
infoinput   = inputdlg({'Measurement Number ='}, ':)', [1 50], {'0'});
nMeasure    = cell2mat(infoinput);


MonitorName    = 'OmniStudioX31-2_Gen2_NoCal_01_40cm';
% MonitorName    = 'OmniStudioX31-2_Gen2_NoCal2';
% MonitorName    = 'OmniStudioX31-2_Gen1_NoCal';
% MonitorName    = 'OmniStudioX31-2_Gen2_NoCal';
% MonitorName    = 'SurfacePro8_NoCal';
OutFileName    = strcat(MonitorName, '_Spectrum', nMeasure, '.txt');
ResultFileName = strcat(MonitorName, '_Result',   nMeasure, '.csv');
xlsxFileName   = strrep(OutFileName, '.txt', '.xlsx');

%% --- Wavelength setup ---
nWaves      = 81;
wavelengths = linspace(380, 380 + 5*(nWaves-1), nWaves);

%% --- Initialize PR-670 ---
CMCheckInit(5);

%% --- Prepare output matrix ---
nGuns      = 3;
nCols      = 1 + 2*nGuns;
outputData = zeros(nWaves, nCols);
outputData(:,1) = wavelengths';

%% --- Pack all shared state into one struct ---
state.gunNames     = {'Red', 'Green', 'Blue'};
state.swatchColors = {[0.8 0.1 0.1], [0.1 0.7 0.1], [0.1 0.3 0.9]};
state.currentGun   = 1;
state.nGuns        = nGuns;
state.nWaves       = nWaves;
state.outputData   = outputData;

%% --- Build GUI ---
fig = uifigure('Name', 'PR-670 RGB Measurement', ...
               'Position', [400 300 480 340], ...
               'Color', [0.12 0.12 0.12], ...
               'Resize', 'off');

uilabel(fig, ...
    'Text', 'PR-670 Primary Measurement', ...
    'Position', [20 290 440 35], ...
    'FontSize', 18, 'FontWeight', 'bold', ...
    'FontColor', [1 1 1], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.12 0.12 0.12]);

instrLabel = uilabel(fig, ...
    'Text', '', ...
    'Position', [20 240 440 40], ...
    'FontSize', 13, ...
    'FontColor', [0.85 0.85 0.85], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.12 0.12 0.12], ...
    'WordWrap', 'on');

swatchPanel = uipanel(fig, ...
    'Position', [190 170 100 50], ...
    'BackgroundColor', [0.3 0.3 0.3], ...
    'BorderType', 'none');

statusLabel = uilabel(fig, ...
    'Text', 'Ready', ...
    'Position', [20 130 440 30], ...
    'FontSize', 12, ...
    'FontColor', [0.6 0.9 0.6], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.12 0.12 0.12]);

progressLabel = uilabel(fig, ...
    'Text', 'Step 0 of 3', ...
    'Position', [20 100 440 25], ...
    'FontSize', 11, ...
    'FontColor', [0.5 0.5 0.5], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.12 0.12 0.12]);

measureBtn = uibutton(fig, ...
    'Text', 'Take Measurement', ...
    'Position', [140 45 200 42], ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.2 0.5 0.9], ...
    'FontColor', [1 1 1]);

% Store handles in fig so callbacks can reach them
fig.UserData.instrLabel    = instrLabel;
fig.UserData.swatchPanel   = swatchPanel;
fig.UserData.statusLabel   = statusLabel;
fig.UserData.progressLabel = progressLabel;
fig.UserData.measureBtn    = measureBtn;
fig.UserData.state         = state;

% Assign callbacks after UserData is set
measureBtn.ButtonPushedFcn = @(~,~) takeMeasurement(fig);

% Set initial GUI state
updateGUI(fig);

%% --- Wait for all measurements ---
uiwait(fig);
pause(1);

% Retrieve outputData from fig before closing
finalState = fig.UserData.state;
outputData = finalState.outputData;
close(fig);

%% --- Save to text file ---
fileId  = fopen(OutFileName, 'w');
tabsOut = repmat('%f\t', 1, nCols-1);
for row = 1:nWaves
    fprintf(fileId, [tabsOut, '%f\n'], outputData(row,:));
end
fclose(fileId);
fprintf('Text data saved to: %s\n', OutFileName);

%% --- Save to xlsx: [wavelength | R_spd | G_spd | B_spd] ---
xlsxData = [outputData(:,1), outputData(:,2), outputData(:,4), outputData(:,6)];
writematrix(xlsxData, xlsxFileName);
fprintf('XLSX saved to: %s\n', xlsxFileName);

%% --- Plot primaries ---
figure;
plot(xlsxData(:,1), xlsxData(:,2), 'r', 'LineWidth', 2); hold on;
plot(xlsxData(:,1), xlsxData(:,3), 'g', 'LineWidth', 2);
plot(xlsxData(:,1), xlsxData(:,4), 'b', 'LineWidth', 2);
xlabel('Wavelength (nm)'); ylabel('Spectral Power');
title('Display Primaries — PR-670'); legend('R','G','B'); grid on; hold off;


%% --- LMS2RGB ---
load SSconefund_0.mat; conefund = SSconefund_0;
convDirec = 2;
spectrum = xlsread(xlsxFileName);

inputDirec1 = [1 0 0];
inputDirec2 = [0 1 0];
inputDirec3 = [0 0 1];

[outDirec1, VecProp1] = NineNumberCalc(conefund, inputDirec1, convDirec, spectrum);
[outDirec2, VecProp2] = NineNumberCalc(conefund, inputDirec2, convDirec, spectrum);
[outDirec3, VecProp3] = NineNumberCalc(conefund, inputDirec3, convDirec, spectrum);

%% --- Save results table ---
L       = [inputDirec1(1); inputDirec2(1); inputDirec3(1)];
M       = [inputDirec1(2); inputDirec2(2); inputDirec3(2)];
S       = [inputDirec1(3); inputDirec2(3); inputDirec3(3)];
VecProp = [VecProp1; VecProp2; VecProp3];

R_old = 2.*[outDirec1(1); outDirec2(1); outDirec3(1)] - 1;
G_old = 2.*[outDirec1(2); outDirec2(2); outDirec3(2)] - 1;
B_old = 2.*[outDirec1(3); outDirec2(3); outDirec3(3)] - 1;
R_new = [outDirec1(1); outDirec2(1); outDirec3(1)];
G_new = [outDirec1(2); outDirec2(2); outDirec3(2)];
B_new = [outDirec1(3); outDirec2(3); outDirec3(3)];

T = table(L, M, S, VecProp, R_old, G_old, B_old, R_new, G_new, B_new);
writetable(T, ResultFileName);
fprintf('Results saved to: %s\n', ResultFileName);


%% =========================================================
%  LOCAL FUNCTIONS  (must be at bottom of script file)
%% =========================================================

function takeMeasurement(fig)
    s     = fig.UserData.state;
    btn   = fig.UserData.measureBtn;
    stat  = fig.UserData.statusLabel;
    prog  = fig.UserData.progressLabel;
    instr = fig.UserData.instrLabel;

    if s.currentGun > s.nGuns
        return;
    end

    btn.Enable  = 'off';
    stat.Text   = sprintf('Measuring %s... please wait', s.gunNames{s.currentGun});
    stat.FontColor = [0.9 0.85 0.3];
    drawnow;

    Beeper;
    [spd, qual] = PR670measspd;
    qual_col    = repmat(qual, s.nWaves, 1);

    s.outputData(:, 2*s.currentGun)     = spd;
    s.outputData(:, 2*s.currentGun + 1) = qual_col;

    stat.Text      = sprintf('%s done. Quality code: %d', s.gunNames{s.currentGun}, qual);
    stat.FontColor = [0.6 0.9 0.6];

    s.currentGun = s.currentGun + 1;

    % Write state back
    fig.UserData.state = s;

    if s.currentGun <= s.nGuns
        updateGUI(fig);
        btn.Enable = 'on';
    else
        instr.Text         = 'All 3 primaries measured. Saving and analysing...';
        prog.Text          = 'Step 3 of 3 — Complete';
        btn.Text           = 'Done';
        btn.BackgroundColor = [0.2 0.7 0.3];
        drawnow;
        uiresume(fig);
    end
end

function updateGUI(fig)
    s     = fig.UserData.state;
    g     = s.currentGun;
    instr = fig.UserData.instrLabel;
    swatch= fig.UserData.swatchPanel;
    stat  = fig.UserData.statusLabel;
    prog  = fig.UserData.progressLabel;
    btn   = fig.UserData.measureBtn;

    instr.Text              = sprintf('Set your display to full %s\nthen click Take Measurement.', s.gunNames{g});
    swatch.BackgroundColor  = s.swatchColors{g};
    prog.Text               = sprintf('Step %d of 3 — measuring %s', g, s.gunNames{g});
    stat.Text               = 'Ready';
    stat.FontColor          = [0.6 0.9 0.6];
    btn.Enable              = 'on';
end