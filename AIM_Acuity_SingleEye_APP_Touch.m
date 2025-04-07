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
function [output_graph, expDuration, testEye, testSName, angleFitobject, AIMcompared, trialRecord] =...
    AIM_Acuity_SingleEye_APP_Touch(eye, view_distance, Patient_Name, number_of_trials, correction_type)

% clear all; close all; commandwindow; % clear all memory, close windows and type into command window
% eye=1; view_distance='40'; Patient_Name='hvh'; Patient_ID = 1; number_of_trials= 1; correction_type = 0;


fprintf("#######################AIM_Acuity_Single_Eye Starting_%s#######################\n", datestr(now,'mm-dd-yyyy HH-MM-SS'))
examtime= datetime;
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
%% read in variable testing parameters:
prompt = {'Subject Initials', 'Weber Contrast','Background Luminance (0-255)','Ring Luminance (0-255)', 'Ring Thickness (deg)', 'Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (40, 300, 400, 600cm)','With(1) of Without(0) Correction','Which eye? left=1,right=2,both=3'};
dlg_title = 'AIM Acuity';
num_lines = 1;

sName=char(Patient_Name); % assign experimenter values to experimnet
tContrast=str2num(char('0.97')); % Weber contrast of letters
meanLUT=str2num(char('194')); % Weber contrast of letters
ringLineDeg=str2num(char('0.2')); % Weber contrast of letters
cSize=str2num(char('1.3')); % diameterof cell rings
nRows=str2num(char('4')); % $ rows and columns on each chart
nCols=str2num(char('4'));
nTrials=number_of_trials; % # charts to run '2'

% scrnWidthCm=str2num(char('70')); % screen width (cm)
[width, height]=Screen('DisplaySize', whichScreen); % adaptive screen width measurement in mm
scrnWidthCm = width/10; % converting mm to cm

viewDistCm=str2num(char(view_distance)); % observer's distance from screen
correction=correction_type; %with (1) or without (2) spectacle correction
testEye = eye; % eye sleected by user 1=left, 2=right,  3=both
ranErrorDeg=90; % random error variable for AIM acuity
cellSizeDeg=5;% hard coded for Nicos experiment, force it to make the right cell size regardless of GUI
nextSizeDeg=4; % size of next chart button

%% save data
if correction ==1
    correction='wCorr';
elseif correction ==0
    correction='woCorr';
end


%which eye tested?
if testEye==1
    testEye='OS';
    txteye='Please cover your right eye and tap to continue';
    soundFile = 'female_cover_right.mp3';
elseif testEye==2
    testEye='OD';
    txteye='Please cover your left eye and tap to continue';
    soundFile = 'female_cover_left.mp3';
elseif testEye==3
    testEye='OU';
    txteye='Look at the screen with both eyes open.Tap to continue';
    soundFile = 'female_use_both_eyes.mp3';
end

testSName=[sName,'_' testEye '_' correction '_' datestr(now,'mm-dd-yyyy HH-MM-SS') ];
dataFile=fullfile(cd,'Patient Data', sprintf('%s_AIM_Acuity.mat', testSName));  % file path to unique Mat files

%% screen and stimulus parameters
rng('default'); % seed random number generator
gammaVal=2.2; % gamma on this system
LMean=255*(meanLUT./255).^(1/gammaVal); % Gamma-corrected background
minSRange=0.2; % minimum letter size range equal 2 logMAR lines

ringLUT=[LMean*0.7 LMean*0.7 LMean*0.7]
targLum=round(meanLUT-(tContrast.*meanLUT)); % luminance of target 0-255
targLum=255*(targLum./255).^(1/gammaVal); % Gamma-corrected target
respLum=255*(ringLUT./255).^(1/gammaVal); % Gamma-corrected target

% feedbackArcRGB=[255 255 255; 160 160 200; 15 15 20]; % bluish ring
% feedbackArcRGB=[160 160 200; 160 160 200; 15 15 20]; % bluish ring
feedbackArcRGB=[LMean*0.5 LMean*0.5 LMean*0.5]; % Grayish ring (120 % times LMean)
feedbackArcRGB=round(255*(feedbackArcRGB./255).^(1/gammaVal)); % Gamma-corrected target

gapArcRGB=[LMean*0.1 LMean*0.1 LMean*0.1]; % 50 % of LMean
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
    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,100); % put some instructions on screen
    DrawFormattedText(windowPtr, txteye, 'center', 'center', 0, LMean(1)); % put this message on screen
    Screen('Flip', windowPtr); % flip to the information screen
    instructionCommand = audioread(soundFile); % sound file when mouse is clicked
    sound(instructionCommand, 22000);
    WaitSecs(5);

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
    end

    instruction = 'Tap where you see the gap of each small black ring';
    instructImage = imread("AIM_Acuity_Image_rotate.jpg");
    Texture = Screen('MakeTexture', windowPtr, instructImage);
    rect = [1050, 1100, 1900, 1500];
    Screen('DrawTexture', windowPtr, Texture, [], rect);

    [~, ~, textBox] = DrawFormattedText(windowPtr,  instruction, 'center', 'center',0, LMean(1)); % give instructions to subject
    Screen('Flip', windowPtr); % flip to the information screen
    command = audioread('female_AIMAcuitytaskinstruction.mp3'); % sound file when mouse is clicked
    sound(command, 22000);
    WaitSecs(5);

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

                                        beepData = audioread('pop.wav'); % sound file when mouse is clicked
                                        sound(beepData, 48000);

                                        Screen('FrameArc',windowPtr,feedbackArcRGB(contNum,:),destRect,mouseOriWRTChosenBox+testGapWidthDeg/2,360-testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw stimulus reported arc

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
                                    %computing 'next' buttom
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
            %% Were the responses significantly different from random chances? If so continue, otherwise show the following screen
            temp_y_90(1:length(trialRecord(trialNo,contNum).oriErr(:)))=ranErrorDeg;
            [h,p]=  ttest(abs(trialRecord(trialNo,contNum).oriErr(:)),temp_y_90')

            if h  == 0 %  minimum response error is not significantly different from random error either 45 or 90 deg
                % Note: random responses, not different from 90/45
                if(trialNo< nTrials)
                    instcommand1 = audioread('female_AIMAcuity_warning.mp3'); % sound file when mouse is clicked
                    sound(instcommand1, 22000)
                    WaitSecs(7)
                end
                trialRecord(trialNo,contNum).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response

            elseif h  == 1 %minimum response error is significantly different from random error either 45 or 90 deg
                % min error is significantly different but we are not sure
                % if its larger than 45/90 deg therefore we do t test
                [h1,p1]=  ttest(abs(trialRecord(trialNo,contNum).oriErr(:)),ranErrorDeg, 'Tail','right')
                if h1 == 0 %  min error data is not signifantly larger than 90/45 degree based onin right sided tailed t test
                    %Note: regular, valid trial
                    if trialNo==1 && trialNo~=nTrials % if the trial is greater than the 1st trial and smaller than the last trial, it will show the Good job screen
                        % Outputs
                        feedbackTxt = 'Well done! The next chart is getting harder.';
                        %                     Speak(feedbackTxt);
                        instcommand3 = audioread('female_welldone_new.mp3'); % sound file when mouse is clicked
                        sound(instcommand3, 22000)

                        WaitSecs(4)
                    elseif trialNo >= 2 && trialNo < (nTrials-1)
                        feedbackTxt = 'You are doing great! Keep going.';
                        instcommand3 = audioread('female_youaredoinggrest_new.mp3'); % sound file when mouse is clicked
                        sound(instcommand3, 22000)
                        WaitSecs(4)
                        %                     Speak(feedbackTxt);
                    elseif trialNo == (nTrials-1)
                        feedbackTxt = 'One more chart left! Do your best.';
                        %                     Speak(feedbackTxt);
                        instcommand4 = audioread('female_onemorechartleft_new.mp3'); % sound file when mouse is clicked
                        sound(instcommand4, 22000)
                        WaitSecs(4)
                    end
                    trialRecord(trialNo,contNum).randomchartresp=0; %count number of random responses for each chart %0 == valid response
                elseif h1 == 1  %  min error data is  signifantly larger than 90/45 degree based onin right sided tailed t test
                    %Note: invalid trial
                    if(trialNo< nTrials)
                        instcommand1 = audioread('female_AIMAcuity_warning.mp3'); % sound file when mouse is clicked
                        sound(instcommand1, 22000)
                        WaitSecs(7)
                    end
                    trialRecord(trialNo,contNum).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response
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
    [angleFitobject,~,angThreshAcuity,AIMci,~,AIMcompared]=fitPFuncOriErr_Aim_Acuity(testLevel',abs(respErr)',90,1,testEye,sName); % fit cumulative Gaussian to the ori error data

    
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

    %     testSName=[sName,'_' testEye '_' correction '_'  datestr(now,'mm-dd-yyyy HH-MM-SS')];

    dataFile=fullfile(cd,'Patient Data', sprintf('%s_AIM_AcuityFailedRun.mat', testSName));  % file path to unique Mat files


    %     save(dataFile,'trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR');
    save(dataFile,'examtime','trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR');
end
end

function [ fitobject, gof, thresholdEst, ci, chiSqFitReject, AIMcompared] = fitPFuncOriErr_Aim_Acuity( tOri, mErr, ranErr, graphData,eye,sName, criterion,AIMtype  )
% Fit data from angular error (target - indicated responses) using cumulative gaussian psychometric function.
%% Input
% tOri= tested orientation angles
% mErr= measured error angle
% ranErr= max error angle ( C=90 in each direction and gratings/bipartie  section =45 )
% graphData= plot figure with outcomes
% eye= which eye has been tested? OS, OD, OU
% sName= subject name
% criterion= which point at the function should be choosen? depending on criterion e.g. 60% correct== 72 degree error in each direction
% AIMtype= function will be fitted without and with constraints, which
% may vary between AIM tests, AIMtype will indicated which fixed variables should be choosen
%% Output
% fitobject= Function returns fitted parameters for semi and constraint fits
% gof= goodness of fit for semi and constraint fits
% thresholdEst=threshold estimates for semi and constraint fits
% ci=95 % Confidence intervals for the fitted parameters for semi and constraint fits
% chiSqFitReject = to test whether data are indeed signifcantly different from zero for semi and constraint fits
% AIMcompared= Contains AIM semi constraint and its comparison e.g. for visual acuity AIM VA and ETDRS equivalent (fully constraint)

%% default values
if nargin<4
    graphData=1;%default plot ON (1) or OFF(~1)
end
if nargin<5
    eye='';%default eye
end
if nargin<6
    sName='XYZ';%default name
end
if nargin<7
    criterion=72;%default etdrs error in degree (40%, 20% in each direction)
end

if nargin<8
    AIMtype=0;%default (0) uses contraint for AIM Acuity, i.e. slopes are 0.13 and noise is 15 (degree)
end
%% AIM type fixed parameters
if AIMtype== 0 %% 0==AIM Acuity, slope= 0.13 and min error=15
    slope_fixed= 0.13;
    noise_fixed= 15;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) AIM parameters: fit semi-constraint function  i.e. boundaries constraint
UBslopesemiconstraint=0.5; % was range(tori)i.e. largest to smallest, but this caused outlier values, hence 0.5 slope which is 3.1fold change over one step
ft = fittype( @(thresh, slope, minErr, x)(minErr + (ranErr-minErr).*(0.5-0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
[fitobject.semiconstraint, gof.semiconstraint] = fit(tOri, mErr,ft, 'StartPoint', [mean(tOri) std(tOri) min(mErr)],  'Upper',[max(tOri) UBslopesemiconstraint ranErr*0.5],  'Lower',[min(tOri) min(diff(tOri)) 1]);
thresholdEst.semiconstraint=fitobject.semiconstraint.thresh;
ci.semiconstraint = confint(fitobject.semiconstraint);
AIMcompared.semiconstraint.logMAR=(sqrt(2).*fitobject.semiconstraint.slope)*(erfinv(((criterion-fitobject.semiconstraint.minErr )/((90-fitobject.semiconstraint.minErr))-0.5)/-0.5))+thresholdEst.semiconstraint; % ang error func
expVals=fitobject.semiconstraint(tOri); % evaluate fit at test levels
chiSqFitReject.semiconstraint=chi2gof(sum(((expVals-mErr).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tOri), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)


%% 2) fit fully constraint function with 1 free parameters, i.e. threshold while slope and min error are fixed
UBslopeconstraint=0.5; % was range(tori)i.e. largest to smallest, but this caused outlier values, hence 0.5 slope which is 3.1fold change over one step
ft = fittype( @(thresh, slope, minErr, x)(minErr + (ranErr-minErr).*(0.5-0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
[fitobject.constraint, gof.constraint] = fit(tOri, mErr,ft, 'StartPoint', [mean(tOri) slope_fixed noise_fixed],  'Upper',[max(tOri) slope_fixed noise_fixed],  'Lower',[min(tOri) slope_fixed noise_fixed]);
thresholdEst.constraint=fitobject.constraint.thresh;
ci.constraint = confint(fitobject.constraint);
AIMcompared.constraint.logMAR=(sqrt(2).*fitobject.constraint.slope)*(erfinv(((criterion-fitobject.constraint.minErr )/((90-fitobject.constraint.minErr))-0.5)/-0.5))+thresholdEst.constraint; % ang error func
expVals=fitobject.constraint(tOri); % evaluate fit at test levels
chiSqFitReject.constraint=chi2gof(sum(((expVals-mErr).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tOri), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphData==1
    %figure()
    %% plot AIM semiconstraint fit
    tOriSort=sort(tOri);
    ci95 = predint(fitobject.semiconstraint,tOriSort, 0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    scatter(tOri, mErr);
    hold on
    ylim([0,ranErr*2]);
    plot(fitobject.semiconstraint, '-');
    h.LineWidth=1;
    h.Color=[0.9290, 0.6940, 0.1250];
    if ~any(isnan(ci95(:))) % no NaNs in confidence intervals
        plot(tOriSort,ci95','--', 'Color',[0, 0.4470, 0.7410],'LineWidth',0.25); %  95% CIs
    end
    hold off
    legend HIDE
    box off
    %     set(gca, 'XScale', 'log');
    ylabel('Angular Error [\circ]','FontSize',12); % label y axis
    ylim([0 180]); % fix upper and lower bounds
    if AIMtype== 0 % 0==AIM Acuity
        %% AIM Acuity's own value is 60% correct with semi-constraints
        AIMcompared.semiconstraint.Feet=20*(10^AIMcompared.semiconstraint.logMAR); % logMAR to 20/X for AIM Acuity VA
        AIMcompared.semiconstraint.Decimal=20/(AIMcompared.semiconstraint.Feet); % Feet to Decimal for AIM Acuity VA
        AIMcompared.semiconstraint.Meter=0.3048*(AIMcompared.semiconstraint.Feet); % Feet to Meters for AIM Acuity VA
        %compare to fully constraint ETDRS equivalent

        AIMcompared.constraint.Feet=20*(10^AIMcompared.constraint.logMAR); % logMAR to 20/X for ETDRS equivalent(72 degree error for C orientation)
        AIMcompared.constraint.Decimal=20/(AIMcompared.constraint.Feet); % Feet to Decimal for AIM Acuity VA
        AIMcompared.constraint.Meter=0.3048*(AIMcompared.constraint.Feet); % Feet to Meters for AIM Acuity VA
        title([sprintf('%s', eye)])
        %         title([sprintf('  %s| ID:%s | AIM VA=  %1.2f logMAR|  20/%.1f | Equi_E_T_D_R_S=  %1.2f logMAR|  20/%.1f', eye,sName,AIMcompared.semiconstraint.logMAR, AIMcompared.semiconstraint.Feet,AIMcompared.constraint.logMAR, AIMcompared.constraint.Feet)],'FontSize',10);
        %         subtitle([' Noise= ' num2str(fitobject.semiconstraint.minErr,2), '\circ ','| Deteriation= ',num2str(fitobject.semiconstraint.slope,2)]);
        xlabel('Letter Size [logMAR]','FontSize',18); % label the x axis
%         xlim([-0.5 1.5]); % fix upper and lower bounds
    else % if nothing indicated
        xlabel('Test Level','FontSize',12); % label the x axis
    end
end
end


