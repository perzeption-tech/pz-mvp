% AIM Acuity
% Program to measure visual acuity based on orientation identification of rotated Landolt C stimulus
% Presents sequence of grids containing random orientation C stimuli of a  range of sizes
% Size of stimuli based on threshold estimate -2 std dev to + 4 stdevs
% Observer clicks perceived orientation of each C and can change response or leave cells empty - scored as incorrect with random orientation error (90 deg)
% Responses classified as 8AFC and fit with cumulative gaussian, mean and stnadard error defines size range for nxt chart
% Data also fit with error function
% Update: Glare ring included, and luminance made more variable
% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021

%% Version 28/01/2022
% function [output_graph, expDuration, testEye, testSName, angleFitobject, AIMcompared] =...
%     AIM_Acuity_SingleEye_APP_Touch(eye, view_distance, Patient_Name, Patient_ID, number_of_trials, correction_type)

clear all; close all; commandwindow; % clear all memory, close windows and type into command window
eye=1; view_distance='40'; Patient_Name='hvh'; Patient_ID = 1; number_of_trials= '1'; correction_type = 0;
stop = 0;
save('stop.mat', 'stop')

examtime= datetime ;
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
%% read in variable testing parameters:
prompt = {'Subject Initials', 'Weber Contrast','Background Luminance (0-255)','Ring Luminance (0-255)', 'Ring Thickness (deg)', 'Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (40, 300, 400, 600cm)','With(1) of Without(0) Correction','Which eye? left=1,right=2,both=3'};
dlg_title = 'AIM Acuity';
num_lines = 1;

sName=char(Patient_Name); % assign experimenter values to experimnet
tContrast=str2num(char('0.97')); % Weber contrast of letters
meanLUT=str2num(char('194')); % Weber contrast of letters
ringLUT=str2num(char('100')); % Weber contrast of letters
ringLineDeg=str2num(char('0.2')); % Weber contrast of letters
cSize=str2num(char('1.3')); % diameterof cell rings
nRows=str2num(char('4')); % $ rows and columns on each chart
nCols=str2num(char('4'));
nTrials=str2num(char(number_of_trials)); % # charts to run '2'

% scrnWidthCm=str2num(char('70')); % screen width (cm)
[width, height]=Screen('DisplaySize', whichScreen); % adaptive screen width measurement in mm
scrnWidthCm = width/10; % converting mm to cm

viewDistCm=str2num(char(view_distance)); % observer's distance from screen

correction=str2num(char(correction_type)); %with (1) or without (2) spectacle correction
testEye = eye; % eye sleected by user 1=left, 2=right,  3=both

%% save data
if correction ==1
    correction='wCorr';
elseif correction ==0
    correction='woCorr';
end
if testEye==1
    testEye='OS';
    txteye='Please cover your right eye and click to continue';
elseif testEye==2
    testEye='OD';
    txteye='Please cover your left eye and click to continue';
elseif testEye==3
    testEye='OU';
    txteye='Please use both eyes and click to continue';
end

testSName=[sName,'_' testEye '_' correction '_' datestr(now,'mm-dd-yyyy HH-MM-SS') ];

dataFile=fullfile(cd,'Patient Data', sprintf('%s_AIM_Acuity.mat', testSName));  % file path to unique Mat files
%% change which monitor, size of surrounding ring and 'next' bottom to use
if viewDistCm<100 % use small HD screen
    disp('IMPORTANT:Use only near AIM acuity < 100cm with additional screen!')
    whichScreen=max(Screen('Screens')); % use highest display # for stimuli
    cellSizeDeg=5;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
    nextSizeDeg=4; % size of next chart button
    scrnWidthCm=31.1; %width of small screen
elseif viewDistCm==300 % if the test is @300cm
    if whichScreen ==0
        whichScreen=max(Screen('Screens')); % use highest display # for stimuli
        cellSizeDeg=1.8;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
        nextSizeDeg=1.5; % size of next chart button
    else %in case multiple screens are connected
        whichScreen=1; % use highest display # for stimuli
        cellSizeDeg=1.8;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
        nextSizeDeg=1.5; % size of next chart button
    end
elseif viewDistCm==400 % if the test is @400cm
    if whichScreen ==0
        whichScreen=max(Screen('Screens')); % use highest display # for stimuli
        cellSizeDeg=1.3;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
        nextSizeDeg=1; % size of next chart button
    else %in case multiple screens are connected
        whichScreen=1; % use highest display # for stimuli
        cellSizeDeg=1.3;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
        nextSizeDeg=1; % size of next chart button
    end
elseif viewDistCm==600 % if the test is @600cm
    if whichScreen ==0
        whichScreen=max(Screen('Screens')); % use highest display # for stimuli
        cellSizeDeg=0.8;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
        nextSizeDeg=0.5; % size of next chart button
        ringLineDeg=0.075;
    else %in case multiple screens are connected
        whichScreen=1; % use highest display # for stimuli
        cellSizeDeg=0.8;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
        nextSizeDeg=0.5; % size of next chart button
        ringLineDeg=0.075;
    end
end
%% screen and stimulus parameters
rng('default'); % seed random number generator
gammaVal=2.2; % gamma on this system
LMean=255*(meanLUT./255).^(1/gammaVal); % Gamma-corrected background
minSRange=0.2; % minimum letter size range equal 2 logMAR lines

targLum=round(meanLUT-(tContrast.*meanLUT)); % luminance of target 0-255
targLum=255*(targLum./255).^(1/gammaVal); % Gamma-corrected target
respLum=255*(ringLUT./255).^(1/gammaVal); % Gamma-corrected target

% feedbackArcRGB=[255 255 255; 160 160 200; 15 15 20]; % bluish ring
% feedbackArcRGB=[160 160 200; 160 160 200; 15 15 20]; % bluish ring

feedbackArcRGB=[LMean*1.2 LMean*1.2 LMean*1.2]; % Grayish ring (120 % times LMean)
feedbackArcRGB=round(255*(feedbackArcRGB./255).^(1/gammaVal)); % Gamma-corrected target

gapArcRGB=[LMean*0.5 LMean*0.5 LMean*0.5]; % 50 % of LMean
gapArcRGB=round(255*(gapArcRGB./255).^(1/gammaVal)); % Gamma-corrected target

interCellGapDeg=1; % updated inter-cell gap cjamged from 0.1 to 1 01-22-2023 SHG & Jan
%Note: intercellgap 0.1 for distance AIM whereas 1 degree for tablet,JSk, 01/22/23


%% other parameters
nAFC=8; % angle for error scoring
saveImage=0; % save stimuli or not

%% Initiating Touch Pad device information gathering
dev = min(GetTouchDeviceIndices([], 1));
if ~isempty(dev)
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
end

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    Screen('Preference','VisualDebugLevel', 0);  % Remove pysch toolbox warning (Red screen)
    %     PsychDebugWindowConfiguration(0, 0.5);

    if ismac
        rval = kPsychNeedRetinaResolution;
        [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean(1),[],[],[],[],[],rval);%,[0 0 3839, 2160]   ,[0 0 2800 1800]);
    else
        [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean(1));%,[0 0 3839, 2160]   ,[0 0 2800 1800]);
    end
    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    imSize=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    srcRect=[0 0 imSize imSize]; % rect for this sized cell
    responseRingPenWidth=pixPerDeg*ringLineDeg; % line width for rresponse ring - TODO justify 0.25 deg
    responseRingRect=[0 0 imSize-responseRingPenWidth/2 imSize-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);

    pixSizeArcmin=60/pixPerDeg; % size of eacdh pixel in arcmin
    minLogMAR=log10(pixSizeArcmin); % minimum target size logMAR
    maxLogMAR=log10((pixSizeArcmin*imSize)/12); % maximum target size logMAR
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % 'next' bottom, halfway down the screen, away from the left

    if ~isempty(dev)
        TouchQueueCreate(windowPtr, dev);
        TouchQueueStart(dev);
        KbReleaseWait;
        HideCursor(windowPtr);
    else
        ShowCursor('Arrow'); % show arrow response cursor
    end
%     instruction_1 = 'Click to start distance tracking for your test';
    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,100); % put some instructions on screen
    DrawFormattedText(windowPtr, txteye, 'center', 'center', 0, LMean(1)); % put this message on screen
    WaitSecs(1);
    Screen('Flip', windowPtr); % flip to the information screen
    Speak(txteye);


    if ~isempty(dev)
        evt = TouchEventGet(dev, windowPtr);
        while isempty(evt)
            evt = TouchEventGet(dev, windowPtr);
            if ~isempty(evt)
                if ismember(evt.Type, [1,2,3,4])
                    continue;
                end
            end
        end

    else
        [mx,my,buttons] = GetMouse; % wait for mouse button release before starting, experiment
        while any(buttons) % if already down, wait for release
            [mx,my,buttons] = GetMouse(windowPtr);
        end
        while ~any(buttons) % wait for new press
            [mx,my,buttons] = GetMouse();
        end
    end
    instruction = 'Click on the gap orientation of each C';

    [~, ~, textBox] = DrawFormattedText(windowPtr,  instruction, 'center', 'center',0, LMean(1)); % give instructions to subject
    Screen('Flip', windowPtr); % flip to the information screen
    Speak(instruction);
    WaitSecs(1);

    Screen('TextSize',windowPtr,48);

    nContrasts=length(tContrast); % how many letter contrasts to test
    slopeEsts=2*ones(size(tContrast)); % start with a rough estimate of slope
    minErrEsts=2*ones(size(tContrast)); % start with a rough estimate of internal error
    startAcuityBnds=[-0.3 1]'; % start with logMAR chart range

    nCells=nRows*nCols;
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
    trialSeed=zeros(nTrials, nContrasts);

    for trialNo=1:nTrials % set up data structure for each chart
        for contNum=1:nContrasts % and for each contrast
            trialRecord(trialNo,contNum) = struct('trialSeed',0,'targSizePerLoc',[],'targLocs',[],'targSize', [],'targOri', [],'matchOri', [],'oriErr', [],'respCorrect', [],'stimSeen', [],'chartTime', []); % start empty array
        end
    end

    %% run experiment
    testGapWidthDeg=360/(5*pi); % gap is always 22.9 deg (1/5 of the line width)
    t=1;%counter for speech command
    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials
        my_coordinate=[]; mx_coordinate=[]; buttons_coordinate=[];ci_contNum=[];
        clear mx; clear my;
        counter_for_distance =0;
        for contNum=1:nContrasts % work though the list of test contrasts
            chartStart=tic; % time for each chart
            if trialNo>1 % there has been at least 1 chart, fit all data for next chart
                testLevel=[]; % set blank list of all levels for this contrast
                respCorrect=[]; % blank list of all stimuli seen for this contrast
                for trialSoFar=1:nTrials % work through all trials
                    testLevel=[testLevel trialRecord(trialSoFar,contNum).targSizePerLoc]; % concatenate all test sizes
                    %                     angErr=[angErr trialRecord(trialNo,contNum).oriErr]; % concatenate all orientation errors
                    respCorrect=[respCorrect  trialRecord(trialSoFar,contNum).respCorrect]; % concatenate all responses
                end
                %                 fitObj=fitPFuncOriErr(tLevel,abs(angErr),1); % fit with current orientation errors
                fitObj=fitPFuncnAFC_FInDAcuity( testLevel, respCorrect, nAFC, 0); % use Matlab's fit function to fit cumulative Gaussian to the data
                %                 logThreshEst(contNum)=fitObj.thresh; % store current threshold estimate
                acuityBnds = confint(fitObj, 0.99); % 99% confidence intervals on threshold estimate
                if any(isnan(acuityBnds)) % fit failed to estimate acuity and 95% CIs
                    acuityBnds(1,1)=fitObj.thresh-2*fitObj.slope; % lower bound is threshold - 2 std
                    acuityBnds(2,1)=fitObj.thresh+2*fitObj.slope; % upper bound is threshold + 2 std
                end
                if any(isnan(acuityBnds)) % fit failed to estimate slope
                    acuityBnds=startAcuityBnds; % use defaults
                end
            else
                acuityBnds=startAcuityBnds; % first trial: use defaults
            end

            minTestLevel=max([minLogMAR acuityBnds(1,1)]); % lower 99% CI in arcmin and minimum size 1 pixel
            maxTestLevel=min([maxLogMAR acuityBnds(2,1)]); % uppper 99% CI in arcmin and maximum size within target circle
            sRange=maxTestLevel-minTestLevel;
            if sRange<minSRange % less than 0.1 log range, increase accordingly
                minTestLevel=minTestLevel-(minSRange-sRange)/2;
                maxTestLevel=maxTestLevel+(minSRange-sRange)/2;
            end
            Screen('FillRect', windowPtr, LMean(contNum));
            Screen('Flip', windowPtr); % clear screen

            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(trialNo, contNum).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(trialNo, contNum).trialSeed); % seed the random number generator

            trialRecord(trialNo,contNum).targLocs=Shuffle(1:nCells); % pick random target locations for target sizes
            trialRecord(trialNo,contNum).targSize=linspace(minTestLevel, maxTestLevel, nCells); % log spaced range of sizes between upper and lower test sizes in arcmin
            trialRecord(trialNo,contNum).targSizePerLoc(trialRecord(trialNo,contNum).targLocs)=trialRecord(trialNo,contNum).targSize; % fill target locations with corresponding contrast

            trialRecord(trialNo,contNum).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,contNum).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,contNum).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,contNum).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

            xStimCenter=screenCenter(1)+xFactor*(imSize+interCellGapPix); % x center of each cell on threen for this size target
            yStimCenter=screenCenter(2)+yFactor*(imSize+interCellGapPix); % y center of each cell on threen for this size target

            for cellNum=1:nCells % draw all random orientation rings
                gapWidthPix=(10.^trialRecord(trialNo,contNum).targSizePerLoc(cellNum))/pixSizeArcmin; % convert logMAR to arcmin, then to pixels in this system
                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord(trialNo,contNum).targOri(cellNum)=testAngle; % update record for test angle
                destRect=CenterRectOnPoint([0 0 5*gapWidthPix 5*gapWidthPix ], xStimCenter(cellNum), yStimCenter(cellNum)); % rect for target size in this cell
                Screen('FrameArc',windowPtr,targLum(contNum),destRect,testAngle+testGapWidthDeg/2,360-testGapWidthDeg,gapWidthPix,gapWidthPix); % arc of required size with required gap
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('FrameOval', windowPtr,respLum(contNum),destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
            end


            Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer

            if ~isempty(dev)
                HideCursor(windowPtr);
                TouchEventFlush(dev);
                blobcol = {};
                blobmin = inf;
            else
                ShowCursor('Arrow'); % show arrow response cursor
            end
            clickedExit=0; % reset clicked finished

            while clickedExit==0 % keep checking until subject click exit location

                if isempty(dev)
                    [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response

                    my_coordinate=[my_coordinate; my];
                    mx_coordinate=[mx_coordinate; mx];
                    buttons_coordinate=[buttons_coordinate; buttons];

                    while any(buttons) % if already down, wait for release
                        [mx,my,buttons] = GetMouse(windowPtr);
                    end
                    while ~any(buttons) % wait for new press
                        [mx,my,buttons] = GetMouse(windowPtr);
                    end
                else
                    while TouchEventAvail(dev)
                        evt = TouchEventGet(dev, windowPtr);
                        id = evt.Keycode;
                        if isinf(blobmin)
                            blobmin = id - 1;
                        end
                        id = id - blobmin;

                        v = 1;
                        if ismember(evt.Type, [2])
                            if v > 1 && blobcol{id}.x == evt.MappedX &&  blobcol{id}.y == evt.MappedY
                                continue;
                            else
                                blobcol{id}.x = evt.MappedX;
                                blobcol{id}.y = evt.MappedY;
                                v = v + 1;

                                if ~isempty(blobcol{id})
                                    if id ~=1 && blobcol{id}.x == blobcol{id-1}.x && blobcol{id}.y == blobcol{id-1}.y
                                        continue;
                                    end

                                    mx = blobcol{id}.x;
                                    my = blobcol{id}.y;
                                    my_coordinate=[my_coordinate; my];
                                    mx_coordinate=[mx_coordinate; mx];                                    

                                    if mx > nextRect(1) && sum(trialRecord(trialNo,contNum).stimSeen)==nCells% observer clicked finished word, only when sum of indications made equals nCells
                                        mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                                        [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                                        mouseOriWRTChosenBox=90+atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                                        clickedExit=1;
                                    else %feedback
                                        
                                        mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                                        [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                                        mouseOriWRTChosenBox=90+atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                                        destRect=CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback
                                        Screen('FrameArc',windowPtr,gapArcRGB(contNum,:),destRect,mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                                        %                                 disp(mx);
                                        beepData = audioread('pop.wav'); % sound file when mouse is clicked
                                        sound(beepData, 48000);

                                        Screen('FrameArc',windowPtr,feedbackArcRGB(contNum,:),destRect,mouseOriWRTChosenBox+testGapWidthDeg/2,360-testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw stimulus reported arc
                                        Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                                        
                                        if rem(counter_for_distance, 4) == 0 
                                            distance = load('test.mat', 'distance');                                  
                                            disp(distance)
                                        end
                                        counter_for_distance = counter_for_distance + 1;
                                        
                                    end
                                    dataMouse(trialNo,contNum).my=my_coordinate;
                                    dataMouse(trialNo,contNum).mx=mx_coordinate;

                                    if isempty(dev)
                                        dataMouse(trialNo,contNum).buttons=buttons_coordinate;
                                    end

                                    trialRecord(trialNo,contNum).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
                                    trialRecord(trialNo,contNum).oriErr(respNum)=diff(unwrap([trialRecord(trialNo,contNum).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % calculate difference between actual and reported orientation
                                    trialRecord(trialNo,contNum).stimSeen(respNum)=1; % mark this cell as response made

                                    %                                     disp(sum(trialRecord(trialNo,cotNum).stimSeen))
                                    %computing 'next' button
                                    if sum( trialRecord(trialNo,contNum).stimSeen)==nCells && nTrials~=trialNo %change   of next button to green once all cells have been clicked, say NEXT
                                        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                                        Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                                        [~, ~, textbounds, ~]=DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [], [],[],[],[],[],nextRect);

                                        beepData = audioread('good.wav'); % sound file when mouse is clicked
                                        sound(beepData, 48000);

                                        Screen('TextBackgroundColor', windowPtr, LMean(1)); % change text background colour to default screen color
                                        Screen('TextColor', windowPtr ,[10 10 10]);
                                    elseif  sum( trialRecord(trialNo,contNum).stimSeen)==nCells && nTrials==trialNo %change color of next button to green once all cells have been clicked, say END
                                        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                                        Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'END' to green
                                        [~, ~, textbounds, ~]=DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [], [],[],[],[],[],nextRect);

                                        beepData = audioread('good.wav'); % sound file when mouse is clicked
                                        sound(beepData, 48000);

                                        Screen('TextBackgroundColor', windowPtr, LMean(1)); % change text background colour to default screen color
                                        Screen('TextColor', windowPtr ,[10 10 10]);
                                    end
                                    Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                                    if abs(trialRecord(trialNo,contNum).oriErr(respNum))<(360/(nAFC*2)) % if within nAFC error
                                        trialRecord(trialNo,contNum).respCorrect(respNum)=1; % score correct
                                    else
                                        trialRecord(trialNo,contNum).respCorrect(respNum)=0; % score incorrect
                                    end

                                    trialRecord(trialNo,contNum).chartTime=toc(chartStart);
                                else
                                    % Below threshold: Kill the blob:
                                    blobcol{id} = [];
                                end

                            end
                        end
                    end
                end
            end
        end % end SFs loop
        if saveImage==1 % saving demo image?
            myImfileName=sprintf('AIMAcuity%d.jpg', trialNo); % create filename
            myImage=Screen('GetImage', windowPtr); % grab screnshot
            imwrite(myImage,myImfileName); % write screenshot to image file
        end
    end % end trials loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    if ~isempty(dev)
        TouchQueueStop(dev);
        TouchQueueRelease(dev);
        ShowCursor(windowPtr);
    end
    Screen('CloseAll'); % close all windows
    

    save(dataFile,'dataMouse','trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR'); % save parameters

    %% AIM error function fit
    output_graph= figure('visible','off');
    testLevel=[]; % set blank of all levels
    respErr=[]; % blank of all response errors
    stimSeen=[]; % blank  of all stimuli seen
    for trialSoFar=1:nTrials % work through all trials
        testLevel=[testLevel trialRecord(trialSoFar,1).targSizePerLoc]; % concatenate all test sizes
        respErr=[respErr trialRecord(trialSoFar,1).oriErr]; % concatenate all orientation errors
        stimSeen=[stimSeen trialRecord(trialSoFar,1).stimSeen]; % concatenate all responses made
    end
    respErr(stimSeen==0)=90; % assign error on ignored cells to 90 deg
    [angleFitobject,~,angThreshAcuity,AIMci,~,AIMcompared]=fitPFuncOriErr(testLevel',abs(respErr)',90,1,testEye,sName); % fit cumulative Gaussian to the ori error data
    stop = 1;
    save('stop.mat', 'stop')
    save(dataFile,'examtime','AIMcompared','angleFitobject','respErr','dataMouse','trialRecord','expDuration',...
        'tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','angThreshAcuity');

catch exception
    if ~isempty(dev)
        TouchQueueRelease(dev);
    end

    disp(exception.message);
    for i = 1:numel(exception.stack)
        exception.stack(i)
    end
    Screen('CloseAll');
    expDuration=toc(expStart); % stop timer - how long was experiment?

    testSName=[sName,'_' testEye '_' correction '_'  datestr(now,'mm-dd-yyyy HH-MM-SS')];

    dataFile=fullfile(cd,'Patient Data', sprintf('%s_AIM_AcuityFailedRun.mat', testSName));  % file path to unique Mat files

    stop = 1;
    save('stop.mat', 'stop')
    %     save(dataFile,'trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR');
    save(dataFile,'examtime','trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR');
end
% end
