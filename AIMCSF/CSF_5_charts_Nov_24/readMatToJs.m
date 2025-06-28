% data = load('your_file.mat'); % replace with your actual filename

data = load('/MATLAB Drive/01_OU_woCorr_11-01-2024 14-02-01AIM_CSF2D.mat')
% data = load('/MATLAB Drive/01_OU_woCorr_11-01-2024 14-05-10AIM_CSF2D.mat')
% data = load('/MATLAB Drive/02_OU_woCorr_11-12-2024 14-51-45AIM_CSF2D.mat')
% data = load('/MATLAB Drive/02_OU_woCorr_11-12-2024 14-56-43AIM_CSF2D.mat')
% data = load('/MATLAB Drive/03_OU_woCorr_11-13-2024 10-58-20AIM_CSF2D.mat')
% data = load('/MATLAB Drive/03_OU_woCorr_11-13-2024 11-04-41AIM_CSF2D.mat')
% data = load('/MATLAB Drive/04_OU_woCorr_11-14-2024 12-27-51AIM_CSF2D.mat')
% data = load('/MATLAB Drive/04_OU_woCorr_11-14-2024 12-31-29AIM_CSF2D.mat')
% data = load('/MATLAB Drive/005Name_OU_woCorr_11-18-2024 10-47-15AIM_CSF2D.mat')
% data = load('/MATLAB Drive/005Name_OU_woCorr_11-18-2024 10-52-21AIM_CSF2D.mat')
% data = load('/MATLAB Drive/006_OU_woCorr_11-18-2024 11-22-23AIM_CSF2D.mat')
% data = load('/MATLAB Drive/006_OU_woCorr_11-18-2024 11-26-29AIM_CSF2D.mat')
% data = load('/MATLAB Drive/007 Name_OU_woCorr_11-18-2024 14-33-59AIM_CSF2D.mat')
% data = load('/MATLAB Drive/007 Name_OU_woCorr_11-18-2024 14-37-40AIM_CSF2D.mat')

% If the .mat file contains only one variable, get its name
fields = fieldnames(data);
if numel(fields) == 1
    mainStruct = data.(fields{1});
else
    mainStruct = data; % Convert the entire structure
end

% Convert to JSON
jsonStr = jsonencode(mainStruct, 'PrettyPrint', true);

% Save JSON string to file
fid = fopen('output.json', 'w');
if fid == -1
    error('Cannot create JSON file');
end
fwrite(fid, jsonStr, 'char');
fclose(fid);

disp('Conversion to JSON complete.');