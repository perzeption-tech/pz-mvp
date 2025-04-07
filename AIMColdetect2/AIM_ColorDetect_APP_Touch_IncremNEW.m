% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021

%% Version 01/28/2023
% function [output_graph, testSName, testEye, trialRecord, randCond, expDuration, angleFitobject, AIMcompared, ci95, correction] = AIM_ColorDetect_APP_Touch(view_distance,Patient_Name, eye, number_of_trials,correction_type)
clear all; close all; commandwindow;n=1; % close windows and type into command window
prompt = {'Subject ID', 'With(1) of Without(0) Correction','Which eye? left=1,right=2,both=3'};
dlg_title = 'AIM Cone Contrast';
num_lines = 1;
def = {'XX', '1', '1'};  % some defaul values; 100cm PR comparison hence default
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
Patient_ID=char(answer(1,1)); % assign experimenter values to experimnet
correction_type=str2num(char(answer(2,1)));  %with (1) or without (2) spectacle correction
eye=str2num(char(answer(3,1))); %1=left,2=right, 3=both
TRIAL_DURATION = 25; % 25 seconds


number_of_trials=1;
    view_distance='40';
    %% 

fprintf("#######################AIM_Color_Detect Starting_%s#######################\n", datestr(now,'mm-dd-yyyy HH-MM-SS'))
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
examtime= datetime ;%date and time variable

% sName=char(Patient_Name); % assign experimenter values to experiment
pedContrast=str2num(char('0.2')); % noise
targSizeDeg=str2num(char('2')); % Size of C- 2 degree as in anamoloscope due to cone receptor density profile
cellSizeDeg=str2num(char('5')); % diameter of cell rings, with sufficent space of noise between stimulus
nRows=str2num(char('4')); %  rows  on each chart
nCols=str2num(char('6')); %columns
nTrials=number_of_trials; % # charts to run
ringLineDeg=str2num(char('0.2')); % indication ring width in degree

[width, height]=Screen('DisplaySize', 0);
scrnWidthCm = width/10;  % # screen width converted to (cm)

viewDistCm=str2num(char(view_distance)); % observer's distance from screen,40cm default
correction=correction_type; %with (1) or without (2) spectacle correction
testEye=eye; %1=left,2=right, 3=both

meanLUM=0.5;
rng('default'); % seed random number generator
gammaVal=2.22; % gamma on this system
LMean=round(255*meanLUM.^(1/gammaVal)); % Gamma-corrected background
noiseRateHz=14; % update rate of noise pedestal

nAFC=8; % angle for error scoring
saveImage=0; % save stimuli or not
interCellGapDeg=1;% change from 0.5 to 1 by Jsk, 01/27/2023
ranErrorDeg=90; % random error variable for AIM acuity

% Vector Proportions for each cone threshold level
VecProp_L=0.191539; % for the Surface tablet
VecProp_M=0.228252; % for the Surface tablet
VecProp_S=0.857758; % for the Surface tablet

% VecProp = [VecProp_L, VecProp_M, VecProp_S];

rng('shuffle')
%% which monitor used?
% rgbRatios=[1.0, -0.102998, -0.00317653; -1.0, 0.478886, -0.0109738; 0.122854, -0.197033, 1.0]; % LG system
% rgbRatios=[1.0, -0.123771, -0.0030043; -1.0, 0.4229717, -0.0301763; 0.1591991, -0.1874024, 1.0]; %HP System, only 3 isolating directions
% rgbRatios_raw=[1.0, -0.123771, -0.0030043; -1.0, 0.4229717, -0.0301763; 0.1591991, -0.1874024, 1.0; -0.1591991, 0.1874024, -1.0]; % HP System, added S decrement
rgbRatios_raw=[1	1	1];% -1	0.431193200772238	-0.020940936705664; 0.137459511322705	-0.176953446953666	1]; % cone isolating directions, surface pro 8, without S decrement
%rgbRatios_raw=[1.0, -0.123771, -0.0030043; -1.0, 0.4229717, -0.0301763; 0.1591991, -0.1874024, 1.0]; % Surface pro 8 (02/02/2023), without S decrement

randCond = randperm(size(rgbRatios_raw,1));
rgbRatios=rgbRatios_raw(randCond,:);% randomization of rows of color directions; long(1), mid(2), short wave length(3)
nColors=size(rgbRatios_raw,1);% how many colors to test

%which correction type?
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
%filename
testSName=[Patient_ID,'_' testEye '_' correction ];%filename
dataFile=fullfile(sprintf('%s_AIM_ColDetect.mat', testSName));  % file path to unique Mat files
%initiate the touchpad device info
dev = min(GetTouchDeviceIndices([], 1));
if ~isempty(dev)
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
end

%start of catch portion
try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    Screen('Preference','VisualDebugLevel', 0);  % Remove pysch toolbox warning (Red screen)
    [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean);% , [0 0 2800 1800]);

    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % height of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen, X
    screenCenter(2)=winRect(2)+scrnHeightPix/2;% center of screen, Y
    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    responseRingPenWidth=pixPerDeg*ringLineDeg; % line width for response ring
    responseRingRect=[0 0 cellDiamPix-responseRingPenWidth/2 cellDiamPix-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);
    testGapWidthDeg=360/(5*pi); % gap is always 22.9 deg (1/5 of the line width)

    targSizePix=round(targSizeDeg*pixPerDeg); % size of target on this system
    mylandoltC=makeLandoltC(cellDiamPix, targSizePix); % make a Ladolt C of required size in background
    rndWin=makeDisk(cellDiamPix);
    checkSizePix=8;% 4 old version, jsk july 24 changed to 8always 0.25 deg: round(0.25*pixPerDeg); % size of random noise pedestal checks
    nPrecomputedFrames=14; %changed from 3 to 8, Jsk, 01/28/2023
    srcRect=[0 0 cellDiamPix cellDiamPix]; % rect for this sized cell

    nextSizeDeg=4; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left

    responseRingCol=[LMean*0.8 LMean*0.8 LMean*0.8]; % RGB color of response ring
    responseWidthDeg=10;
    slopeEsts=2*ones(size(targSizeDeg)); % start with a rough estimate of slope
    minErrEsts=2*ones(size(targSizeDeg)); % start with a rough estimate of internal error
      mincontbnd=-2;
    maxcontbnd=0;
    startContrastBnds=[mincontbnd maxcontbnd]'; %-3= .1% to 0= 100%
    
    minCRange=0.5; % minium range of contrasts on screen

    if ~isempty(dev)
        TouchQueueCreate(windowPtr, dev);
        TouchQueueStart(dev);
        KbReleaseWait;
        HideCursor(windowPtr);
    else
        ShowCursor('Arrow'); % show arrow response cursor
    end

    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,100);                               % put some instructions on screen
    %     textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate);

    DrawFormattedText(windowPtr, txteye, 'center', 'center', 0, LMean); % put this message on screen
    %     ShowCursor('Arrow'); % show arrow response cursor
    WaitSecs(1); %UX consideration: give participant time to process the task information
    Screen('Flip', windowPtr); % flip to the information screen
    instructionCommand = audioread(soundFile); % sound file when mouse is clicked
    sound(instructionCommand, 22000);
    WaitSecs(5);

    % begin to check whether a touch occured
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
%     instructions = 'Tap where you see the gap of each colored ring \n otherwise make your best guess';
%     instructImage = imread("AIM_ColorDetect_Image_Rotate.jpg");
%     Texture = Screen('MakeTexture', windowPtr, instructImage);
%     rect = [1050, 1100, 1900, 1500];
% 
%     Screen('DrawTexture', windowPtr, Texture, [], rect);
%     [~, ~, textBox] = DrawFormattedText(windowPtr,  instructions, 'center', 'center',0, LMean); % give instructions to subject
%     Screen('Flip', windowPtr); % flip to the information screen
%     command = audioread('female_AIMColDetect_instructions.mp3'); % sound file when mouse is clicked
%     sound(command, 22000);
%     WaitSecs(6);
%     Screen('TextSize',windowPtr,48);

    nCells=nRows*nCols;
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
    trialSeed=zeros(nTrials, nColors);
    xStimCenter=screenCenter(1)+xFactor*(cellDiamPix+interCellGapPix); % x center of each cell on threen for this size target
    yStimCenter=screenCenter(2)+yFactor*(cellDiamPix+interCellGapPix); % y center of each cell on threen for this size target
    destRect=zeros(nCells, 4);

    for cellNum=1:nCells % work through each cell
        destRect(cellNum, :)=CenterRectOnPoint(srcRect, xStimCenter(cellNum), yStimCenter(cellNum));
    end

    for trialNo=1:nTrials % set up data structure for each chart
        for colNum=1:nColors % and for each contrast
            trialRecord(trialNo,colNum) = struct('trialSeed',0,'targContrastPerLoc',[],'targLocs',[],'targContrast', [],'targOri', [],'matchOri', [],'oriErr', [],'respCorrect', [],'stimSeen', [],'chartTime', []); % start empty array
        end
    end
    chart = 0;

    %% run experiment
    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials
        for colNum=1:nColors % work though the list of test contrasts
            chartStart=tic; % time for each chart
            chart = chart + 1;
            if trialNo>1 % there has been at least 1 chart, fit all data for next chart
                angErr=[]; % set blank list of all levels for this contrast
                testLevel=[]; % set blank list of all levels for this contrast
                stimSeen=[]; % blank list of all stimuli seen for this SF
                respCorrect=[]; % blank list of all stimuli seen for this SF
                for trialSoFar=1:nTrials % work through all trials
                    testLevel=[testLevel trialRecord(trialSoFar,colNum).targContrastPerLoc]; % concatenate all test sizes
                    angErr=[angErr trialRecord(trialSoFar,colNum).oriErr]; % concatenate all orientation errors
                    stimSeen=[stimSeen trialRecord(trialSoFar,colNum).stimSeen]; % concatenate all responses made
                    respCorrect=[respCorrect  trialRecord(trialSoFar,colNum).respCorrect]; % concatenate all orientation errors
                end
                angErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 /45 deg
                %                 fitObj=fitPFuncOriErr_AIMColDetect(testLevel',abs(angErr)',1); % fit with current orientation errors
                fitObj=fitPFuncnAFC_AIM_ColorDetect(testLevel, respCorrect, nAFC, 0); % use matlab's fit function to fit cumulative Gaussian to the data
                csBnds = confint(fitObj, 0.95); % 95% confidence intervals on contrast threshold estimate
                if any(isnan(csBnds)) % fit failed to estimate acuity and 99% CIs
                    contrastBound(1)=fitObj.thresh-2*fitObj.slope; % lower bound is threshold - 2 std
                    contrastBound(2)=fitObj.thresh+2*fitObj.slope; % upper bound is threshold + 2 std
                else
                    contrastBound(1)=csBnds(1);
                    contrastBound(2)=csBnds(2);
                end
                if any(isnan(contrastBound)) % fit failed to estimate slope
                    contrastBound=startContrastBnds; % use defaults
                end
            else
                contrastBound=startContrastBnds; % first trial: use defaults
            end
            %save guard for fit to never breach boundaries of physically possible stimulus properties
            minTestLevel=max([-3 contrastBound(1)]); % lower 99% CI contrast, minimum 0.1%
            minTestLevel=min([log10(1-pedContrast)-minCRange minTestLevel]); % prevent minTestLevel from getting too large (>0)
            % when fitted curve is flat,abs(contrastBound) becomes too large, so force it to be log10(1-pedContrast)-minCRange
            maxTestLevel=min([log10(1-pedContrast) contrastBound(2)]); % uppper 99% CI in arcmin and maximum size within target circle
            if maxTestLevel<minTestLevel
                maxTestLevel=log10(1-pedContrast);
                minTestLevel=log10(1-pedContrast)-minCRange;
            end % prevent maxTestLevel from getting
            cRange=maxTestLevel-minTestLevel;
            if cRange<minCRange % less than 0.1 log range, increase accordingly
                midLevel=minTestLevel+(maxTestLevel-minTestLevel)/2;
                minTestLevel=midLevel-cRange/2;
                maxTestLevel=midLevel+cRange/2;
            end

            vbl= Screen('Flip', windowPtr); % clear screen
            buttonPressTime=vbl;
            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(trialNo, colNum).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(trialNo, colNum).trialSeed); % seed the random number generator

            trialRecord(trialNo,colNum).targLocs=Shuffle(1:nCells); % pick random target locations for target contrasts
            trialRecord(trialNo,colNum).targContrast=linspace(minTestLevel, maxTestLevel, nCells); % lin spaced range of contrasts between upper and lower test contrasts
            trialRecord(trialNo,colNum).targContrastPerLoc(trialRecord(trialNo,colNum).targLocs)=trialRecord(trialNo,colNum).targContrast; % fill target locations with corresponding contrast

            trialRecord(trialNo,colNum).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,colNum).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,colNum).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,colNum).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

            for cellNum=1:nCells % draw all random orientation rings
                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord(trialNo,colNum).targOri(cellNum)=testAngle; % update record for test angle
                myTarg=imrotate(mylandoltC, -testAngle,'bilinear', 'crop'); % rotate the C by required amount

                for preFrameNum=1:nPrecomputedFrames % make as manytextures as requested
                    myNoiseRGB=zeros([cellDiamPix cellDiamPix 3]); % blank stimulus
                    myNoise=randn(round(cellDiamPix/checkSizePix)); % random downsampled noise
                    if checkSizePix>1 % need to resize the noise
                        myNoise=imresize(myNoise, [cellDiamPix cellDiamPix], 'nearest');
                    end
                    myNoise=myNoise/max(abs(myNoise(:))); % scale -1 to +1
                    myNoise=pedContrast*(myNoise.*rndWin); % scale noise to pedestal level
                    for rgbLayer=1:3
                        myNoiseRGB(:,:,rgbLayer)=myNoise+myTarg*10.^trialRecord(trialNo,colNum).targContrastPerLoc(cellNum)*rgbRatios(colNum, rgbLayer);
                    end
                    myNoiseRGB=127+127*myNoiseRGB; % scale to mean luminance
                    myNoiseRGB=255*(myNoiseRGB/255).^(1/gammaVal); % gamma correct
                    myNoiseRGB(myNoiseRGB<0)=0; % catch out of range values
                    myNoiseRGB(myNoiseRGB>255)=255;
                    myNoiseFloor=floor(myNoiseRGB); % base LUT value for each pixel
                    myNoiseResidual=myNoiseRGB-myNoiseFloor; % residual floating point
                    MyNoiseBit=binornd(1, myNoiseResidual)+myNoiseFloor; % binomial random sample with probability from residual
                    myTex(cellNum, preFrameNum)=Screen('MakeTexture', windowPtr,MyNoiseBit);%, [],[],2); % write image to texture, scaled
                end
            end

            if ~isempty(dev)
                HideCursor(windowPtr);
                TouchEventFlush(dev);
                blobcol = {};
                blobmin = inf;

            else
                ShowCursor('Arrow'); % show arrow response cursor
            end

            clickedExit=0; % reset clicked finished
            frameNum=0;
            %%%%%%%%%%%%%%%%%%%%

            % Start the trial timer
            trialStartTime = tic;

            clickedExit = 0;
            frameNum = 0;
            [mx, my, buttons] = GetMouse(windowPtr);

            while toc(trialStartTime) < TRIAL_DURATION && clickedExit == 0
                frameNum = frameNum + 1;

                % Draw stimuli
                for cellNum = 1:nCells
                    Screen('DrawTexture', windowPtr, myTex(cellNum,1+mod(frameNum,nPrecomputedFrames)), srcRect, destRect(cellNum,:));
                    Screen('FrameOval', windowPtr, responseRingCol, destRect(cellNum,:), responseRingPenWidth);
                    if trialRecord(trialNo,colNum).stimSeen(cellNum) == 1
                        Screen('FrameOval', windowPtr, [LMean*0.9 LMean*0.9 LMean*0.9], destRect(cellNum,:), responseRingPenWidth);
                        Screen('FrameArc', windowPtr, [LMean*1.2 LMean*1.2 LMean*1.2], destRect(cellNum,:), 90+trialRecord(trialNo,colNum).matchOri(cellNum)-responseWidthDeg/2, responseWidthDeg, responseRingPenWidth, responseRingPenWidth);
                    end
                end

                vbl = Screen('Flip', windowPtr, vbl+(1/noiseRateHz));

                % Check for touch events
                if ~isempty(dev)
                    while TouchEventAvail(dev)
                        evt = TouchEventGet(dev, windowPtr);
                        ismember(evt.Type,[1,2])%touch down or move
                        mx=evt.MappedX;
                        my=evt.MappedY;
                        mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                        [~,respNum]=min(mouseDistFromEachBox(:));

                        if mouseDistFromEachBox(respNum) <= cellDiamPix

                            trialRecord(trialNo,colNum).matchOri(respNum)=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                            trialRecord(trialNo,colNum).oriErr(respNum)=diff(unwrap([trialRecord(trialNo,colNum).targOri(respNum),trialRecord(trialNo,colNum).matchOri(respNum)]/180*pi)*180/pi);% calculate difference betwen actual and reporte orientation
                            %                         trialRecord(trialNo,colNum).oriErrSorted(trialRecord(trialNo,colNum).targRanLocs(respNum))=trialRecord(trialNo,colNum).oriErr(respNum);
                            trialRecord(trialNo,colNum).stimSeen(respNum)=1;
                            if abs(trialRecord(trialNo,colNum).oriErr(respNum))<(360/(nAFC*2)) % if within nAFC error
                                trialRecord(trialNo,colNum).respCorrect(respNum)=1; % score correct
                            else
                                trialRecord(trialNo,colNum).respCorrect(respNum)=0; % score incorrect
                            end
%                             beepData = audioread('pop.wav'); % sound file when mouse is clicked
%                             sound(beepData, 50000);
% WaitSecs(0.2)
                        end


                        % Check if exit condition is met
                        if mx > nextRect(1) && mx < nextRect(3) && my > nextRect(2) && my < nextRect(4) && sum(trialRecord(trialNo).stimSeen) == nCells
                            clickedExit = 1;
                            beepData = audioread('good.wav');
                            sound(beepData, 50000);
                            break;
                        end
                    end
                    %end


                else
                    % Mouse input processing (if no touch device)
                    [mx, my, buttons] = GetMouse(windowPtr);
                    if any(buttons)
                        % Process mouse click
                        % (Your existing mouse click processing code here)
                    end
                end
           
                % Check if trial duration has elapsed
                if toc(trialStartTime) >= TRIAL_DURATION
                    clickedExit = 1;
                end
            end





            %%%%%%%%%%%%%%%%%%%%%


            % while clickedExit==0 % keep checking until subject click exit location
            %     frameNum=frameNum+1;
            %     %                 Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
            %     %                 DrawFormattedText(windowPtr, sprintf('%s', 'next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
            %     for cellNum=1:nCells % work through each cell
            %         Screen('DrawTexture', windowPtr, myTex(cellNum,1+mod(frameNum,nPrecomputedFrames)), srcRect,destRect(cellNum,:)); % draw texture to on screen location
            %         Screen('FrameOval',windowPtr,responseRingCol, destRect(cellNum,:), responseRingPenWidth);  % draw box around cells
            %         if trialRecord(trialNo,colNum).stimSeen(cellNum)==1 % this cell clicked - put a ring on it
            %             Screen('FrameOval',windowPtr,[LMean*0.9 LMean*0.9 LMean*0.9], destRect(cellNum,:), responseRingPenWidth);  % draw ring around cells
            %             Screen('FrameArc',windowPtr,[LMean*1.2 LMean*1.2  LMean*1.2],destRect(cellNum,:),90+trialRecord(trialNo,colNum).matchOri(cellNum)-responseWidthDeg/2,responseWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc - add 90 because PTB 0 = up
            %         end
            %     end
            %     vbl= Screen('Flip', windowPtr, vbl+(1/noiseRateHz)); % show stimulus after 1/?Hz
            %     %                 vbl= Screen('Flip', windowPtr, vbl+0.1); % show stimulus after 200msec
            %
            %     if isempty(dev)
            %         [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
            %     else
            %         while TouchEventAvail(dev)
            %             evt = TouchEventGet(dev, windowPtr);
            %             id = evt.Keycode;
            %             if isinf(blobmin)
            %                 blobmin = id - 1;
            %             end
            %             id = id - blobmin;
            %
            %             v = 1;
            %             if ismember(evt.Type, [2])
            %                 if evt.MappedX < destRect(16, 3) && evt.MappedX > destRect(1,1) && evt.MappedY < destRect(16, 4) && evt.MappedY > destRect(1,2)
            %                     beepData = audioread('pop.wav'); % sound file when mouse is clicked
            %                     sound(beepData, 50000);
            %                 end
            %                 if v > 1 && blobcol{id}.x == evt.MappedX &&  blobcol{id}.y == evt.MappedY
            %                     continue;
            %                 else
            %                     blobcol{id}.x = evt.MappedX;
            %                     blobcol{id}.y = evt.MappedY;
            %                     v = v + 1;
            %                     %% timer to automatically stop the research
            %                     trialStartTime = tic; % Start the trial timer
            %
            %                     clickedExit = 0; % reset clicked finished
            %
            %                     while toc(trialStartTime) < TRIAL_DURATION && ~clickedExit
            %                         [mx, my, buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
            %                         while any(buttons) % if already down, wait for release
            %                             [mx, my, buttons] = GetMouse(windowPtr);
            %                         end
            %                         while ~any(buttons) % wait for new press
            %                             [mx, my, buttons] = GetMouse(windowPtr);
            %                             if toc(trialStartTime) >= TRIAL_DURATION
            %                                 clickedExit = 1;
            %                                 break;
            %                             end
            %                         end
            %
            %                         if clickedExit
            %                             break;
            %                         end
            %                         if ~isempty(blobcol{id})
            %                             if id ~=1 && blobcol{id}.x == blobcol{id-1}.x && blobcol{id}.y == blobcol{id-1}.y
            %                                 continue;
            %                             end
            %
            %                             mx = blobcol{id}.x;
            %                             my = blobcol{id}.y;
            %
            %                             %                                     [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
            %                             mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
            %                             [~,respNum]=min(mouseDistFromEachBox(:));
            %
            %                             buttonPressTime=vbl; % note screen time when button was pressed
            %
            %
            %                             if mx > nextRect(1) && mx < nextRect(3) && my >nextRect(2) && my< nextRect(4) && sum(trialRecord(trialNo).stimSeen)==nCells % observer clicked finished word and all cells were responded to
            %                                 clickedExit=1;
            %                                 beepData = audioread('good.wav'); % sound file when mouse is clicked
            %                                 sound(beepData, 50000);
            %                             else
            %                                 trialRecord(trialNo,colNum).matchOri(respNum)=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
            %                                 trialRecord(trialNo,colNum).oriErr(respNum)=diff(unwrap([trialRecord(trialNo,colNum).targOri(respNum),trialRecord(trialNo,colNum).matchOri(respNum)]/180*pi)*180/pi);% calculate difference betwen actual and reporte orientation
            %                                 %                         trialRecord(trialNo,colNum).oriErrSorted(trialRecord(trialNo,colNum).targRanLocs(respNum))=trialRecord(trialNo,colNum).oriErr(respNum);
            %                                 trialRecord(trialNo,colNum).stimSeen(respNum)=1;
            %                                 if abs(trialRecord(trialNo,colNum).oriErr(respNum))<(360/(nAFC*2)) % if within nAFC error
            %                                     trialRecord(trialNo,colNum).respCorrect(respNum)=1; % score correct
            %                                 else
            %                                     trialRecord(trialNo,colNum).respCorrect(respNum)=0; % score incorrect
            %                                 end
            %
            %                             end
            %
            %                             trialRecord(trialNo,colNum).chartTime=toc(chartStart);
            %                         else
            %                             blobcol{id} = [];
            %
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %computing 'next' buttom gets green
  
            if sum( trialRecord(trialNo,colNum).stimSeen)==nCells && nTrials==trialNo && colNum==nColors%Create END button once all cells have been clicked for all trials and all colors

                Screen('FillOval',windowPtr,[LMean*0.9 LMean*0.9 LMean*0.9], nextRect);  % draw exit button, green
                Screen('TextBackgroundColor', windowPtr, [LMean*0.9 LMean*0.9 LMean*0.9]); % change text background color of 'NEXT' to green
                DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [], [],[],[],[],[],nextRect); %  char(26)
                Screen('TextBackgroundColor', windowPtr, [LMean*0.9 LMean*0.9 LMean*0.9]); % change text background color of 'NEXT' to green

            elseif sum( trialRecord(trialNo,colNum).stimSeen)==nCells %Create NEXT button once all cells have been clicked for a trial
                Screen('FillOval',windowPtr,[LMean*0.9 LMean*0.9 LMean*0.9], nextRect);  % draw NEXT button, green
                Screen('TextBackgroundColor', windowPtr, [LMean*0.9 LMean*0.9 LMean*0.9]); % change text background color of 'NEXT' to green
                DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [], [],[],[],[],[],nextRect); %  char(26)
                Screen('TextBackgroundColor', windowPtr, [LMean*0.9 LMean*0.9 LMean*0.9]); % change text background color of 'NEXT' to green
            end
%              if saveImage ==1% saving demo image?
%         myImfileName=sprintf('AIM_ColorDetect%d.jpg', colNum); % create filename
%         myImage=Screen('GetImage', windowPtr); % grab screnshot
%         imwrite(myImage,myImfileName); % write screenshot to image file
%     end
        end
    end
    %% Were the responses significantly different from random chances? If so continue, otherwise show the following screen
    temp_y_90(1:length(trialRecord(trialNo,colNum).oriErr(:)))=ranErrorDeg;
    [h,p]=  ttest(abs(trialRecord(trialNo,colNum).oriErr(:)),temp_y_90')

    if h==0 %  minimum response error is not significantly different from random error either 45 or 90 deg
        % Note: random responses, not different from 90/45
        if chart < (nTrials * nColors)
            instcommand1 = audioread('female_makesure_AIMColorDetect.mp3'); % sound file when mouse is clicked
            sound(instcommand1, 22000)
            WaitSecs(5)
        end
        trialRecord(trialNo,colNum).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response

    elseif h==1 %minimum response error is significantly different from random error either 45 or 90 deg
        % min error is significantly different but we are not sure
        % if its larger than 45/90 deg therefore we do t test
        [h1,p1]=  ttest(abs(trialRecord(trialNo,colNum).oriErr(:)),ranErrorDeg, 'Tail','right')
        if h1==0 %  min error data is not signifantly larger than 90/45 degree based onin right sided tailed t test
            %Note: regular, valid trial
            if chart < ((nTrials * nColors)-1)
                instcommand3 = audioread('female_youaredoinggrest_new.mp3'); % sound file when mouse is clicked
                sound(instcommand3, 22000)
                WaitSecs(4)
            elseif chart < (nTrials * nColors)
                instcommand4 = audioread('female_onemorechartleft_new.mp3'); % sound file when mouse is clicked
                sound(instcommand4, 22000)
                WaitSecs(4)
            end
            trialRecord(trialNo,colNum).randomchartresp=0; %count number of random responses for each chart %0 == valid response
        elseif h1 == 1  %  min error data is  signifantly larger than 90/45 degree based onin right sided tailed t test
            %Note: invalid trial
            if chart < (nTrials * nColors)
                instcommand1 = audioread('female_makesure_AIMColorDetect.mp3'); % sound file when mouse is clicked
                sound(instcommand1, 22000)
                WaitSecs(5)
            end
            trialRecord(trialNo,colNum).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response
        end
    end

    trialRecord(trialNo,colNum).chartTime=toc(chartStart);
%     if saveImage ==1% saving demo image?
%         myImfileName=sprintf('AIM_ColorDetect%d.jpg', colNum); % create filename
%         myImage=Screen('GetImage', windowPtr); % grab screnshot
%         imwrite(myImage,myImfileName); % write screenshot to image file
%     end
    %         end % end stimulus type loop
    % end % end trials loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    if ~isempty(dev)
        TouchQueueStop(dev);
        TouchQueueRelease(dev);
        ShowCursor(windowPtr);
    end
    Screen('CloseAll'); % close all windows

    % testSName=[sName,'_' testEye '_' correction '_' num2str(n)];
    % while exist([testSName,'AIM_CS.mat'],'file') ~= 0 % check if data file exists with this subject ID
    %     n=n+1;
    %     testSName=[testSName,num2str(n)]; % if so create new name until unique name found
    % end
%     dataFile=sprintf('%sAIM_CS_Scone.mat', testSName);  % file path to unique Mat files
%     save(dataFile,'examtime','correction','testEye','trialRecord','expDuration','meanLUM','LMean','targSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','randCond','rgbRatios_raw','rgbRatios'); % save parameters

    %     output_graph(1)=figure('visible','off');
    figure
    [~,idx]=sort(randCond);

    for colNum=1:nColors % plot figure with fits for interleaved conditions - attempt to fit all data
        % if idx(colNum)==1 %plot long wave
        subplot(1,nColors,colNum);
        testLevel=[]; % set blank list of all levels for this SF
        respErr=[]; % blank list of all stimuli seen for this SF
        stimSeen=[]; % blank list of all stimuli seen for this SF
        for trialSoFar=1:nTrials % work through all trials
            testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
            respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
            stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
        end
        respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
        testLevel=10.^testLevel;% JH added 11/18/2022
        if colNum==1% Converting the threshold in machine units to cone constrat unit
%             testLevel = VecProp_L*testLevel;
%         elseif colNum==2
%             testLevel = VecProp_M*testLevel;
%         elseif colNum==3
%             testLevel = VecProp_S*testLevel;
        end

        [angleFitobject{colNum},~,angThreshCS{colNum},~,~,AIMcompared{colNum}]=fitPFuncOriErr_AIMColDetect(testLevel',abs(respErr)',ranErrorDeg,1,testEye,Patient_ID);
        if colNum==1
%             subtitle('L-cone Threshold')
%             ci95.longwave = confint(angleFitobject{(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
%         elseif colNum==2
%             subtitle('M-cone Threshold')
%             ci95.midwave = confint(angleFitobject{(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
%         elseif colNum==3
            subtitle('Increment Threshold')
            ci95.shortwave = confint(angleFitobject{(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
        end
    end
    %[angleFitobject,~,angThreshCS]=fitPFuncOriErr_AIMColDetect(testLevel',abs(respErr)',ranErrorDeg,1,testEye,sName,60); % 60 degree==67% 2/3 criterion Pfit cumulative Gaussian to the ori error data

    %             xlabel('Linear Contrast'); % label the x axis;JH11/18/2022:was log10(Contrast), otherwise can not use error function
    %             title(sprintf('%s  | %.2fdeg |threshold %.2f%% |CI(%.2f %.2f)', sName, 100*10.^angThreshCS.unconstraint, 100*10.^csBnds(1), 100*10.^csBnds(2)));
    %             ylim([0 180]); % fix upper and lower bounds

    %         elseif idx(colNum)==2 %plot mid wavelength
    %             subplot(1,nColors,2);
    %             testLevel=[]; % set blank list of all levels for this SF
    %             respErr=[]; % blank list of all stimuli seen for this SF
    %             stimSeen=[]; % blank list of all stimuli seen for this SF
    %             for trialSoFar=1:nTrials % work through all trials
    %                 testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
    %                 respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
    %                 stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
    %             end
    %             respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
    %             testLevel=10.^testLevel;% JH added 11/18/2022
    %             testLevel = testLevel * VecProp_M; % Converting the threshold in machine units to cone constrat unit
    %             [angleFitobject{2},~,angThreshCS{2},~,~,AIMcompared{2}]=fitPFuncOriErr_AIMColDetect(testLevel',abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    %
    %             % [angleFitobject,~,angThreshCS]=fitPFuncOriErr_AIMColDetect(testLevel',abs(respErr)',ranErrorDeg,1,testEye,sName,60); % 60 degree==67% 2/3 criterion Pfit cumulative Gaussian to the ori error data
    %             ci95.midwave = confint(angleFitobject{idx(colNum)}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    %             %             xlabel('Linear Contrast'); % label the x axis;JH11/18/2022:was log10(Contrast), otherwise can not use error function
    %             %             title(sprintf('%s  | %.2fdeg |threshold %.2f%% |CI(%.2f %.2f)', sName, 100*10.^angThreshCS.unconstraint, 100*10.^csBnds(1), 100*10.^csBnds(2)));
    %             %             ylim([0 180]); % fix upper and lower bounds
    %             subtitle(['M-cone Threshold'])
    %
    %         elseif idx(colNum)==3%plot short wavelength
    %             subplot(1,nColors,3);
    %
    %             testLevel=[]; % set blank list of all levels for this SF
    %             respErr=[]; % blank list of all stimuli seen for this SF
    %             stimSeen=[]; % blank list of all stimuli seen for this SF
    %             for trialSoFar=1:nTrials % work through all trials
    %                 testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
    %                 respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
    %                 stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
    %             end
    %             respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
    %             testLevel=10.^testLevel;% JH added 11/18/2022
    %             testLevel = testLevel * VecProp_S; % Converting the threshold in machine units to cone constrat unit
    %             [angleFitobject{3},~,angThreshCS{3},~,~,AIMcompared{3}]=fitPFuncOriErr_AIMColDetect(testLevel',abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    %
    %             %[angleFitobject,~,angThreshCS]=fitPFuncOriErr_AIMColDetect(testLevel',abs(respErr)',ranErrorDeg,1,testEye,sName,60); % 60 degree==67% 2/3 criterion Pfit cumulative Gaussian to the ori error data
    %             ci95.shortwave = confint(angleFitobject{idx(colNum)}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    %             %             xlabel('Linear Contrast'); % label the x axis;JH11/18/2022:was log10(Contrast), otherwise can not use error function
    %             %             ylabel('Angular Error [ \circ ]'); % label y axis
    %             %         title(sprintf('%s  | %.1fdeg |threshold %.2f%% |CI(%.2f %.2f)', sName, 100*10.^angThreshCS{idx(colNum)}.unconstraint, 100*10.^csBnds(1), 100*10.^csBnds(2)));
    %             %             ylim([0 180]); % fix upper and lower bounds
    %             subtitle(['S-cone Threshold'])
    %         end
    %     end

    %% summary threshold plot
    %     output_graph(2) = figure('visible','off');
    figure
    for colorNum=1:nColors % process data foreach color
        if colorNum==1 %long
        %     e= errorbar(colorNum,angThreshCS{colorNum}.unconstraint, angThreshCS{colorNum}.unconstraint-ci95.longwave(1,1),ci95.longwave(2,1)-angThreshCS{colorNum}.unconstraint);hold on
        %     e.Marker = 'o';
        %     e.MarkerSize = 7;
        %     e.Color = 'red';
        %     e.CapSize = 15;
        % elseif colorNum==2 %mid
        %     e= errorbar(colorNum,angThreshCS{colorNum}.unconstraint, angThreshCS{colorNum}.unconstraint-ci95.midwave(1,1),ci95.midwave(2,1)-angThreshCS{colorNum}.unconstraint);hold on
        %     e.Marker = 'o';
        %     e.MarkerSize = 7;
        %     e.Color = 'green';
        %     e.CapSize = 15;
        % elseif colorNum==3 %short
            e= errorbar(colorNum,angThreshCS{colorNum}.unconstraint, angThreshCS{colorNum}.unconstraint-ci95.shortwave(1,1),ci95.shortwave(2,1)-angThreshCS{colorNum}.unconstraint);hold on
            e.Marker = 'o';
            e.MarkerSize = 7;
            e.Color = 'blue';
            e.CapSize = 15;
        end
        title('AIM Color Detection')
        subtitle(['Increm Contrast:' num2str(angThreshCS{colorNum}.unconstraint,3) ], 'FontSize',14);%       M-Cone        S-Cone
        xlabel('Stimulus types');
        xlim([0 nColors+1]);
        ylabel('Threshold Contrast');
        ylim([0 0.2]);%based on prelimnary results FInD Color study
        %         legend([e;e;e], {'Long-Wave', 'Mid-Wave', 'Short-Wave'})
    end

    % testSName=[sName,'_' testEye '_' correction '_' num2str(n)];
    while exist([testSName,'AIM_CS_Increm.mat'],'file') ~= 0 % check if data file exists with this subject ID
        n=n+1;
        testSName=[testSName,num2str(n)]; % if so create new name until unique name found
    end
    dataFile=sprintf('%sAIM_CS_Increm.mat', testSName);  % file path to unique Mat files
    save(dataFile,'AIMcompared','ci95', 'examtime','correction','testEye','trialRecord','angleFitobject', 'angThreshCS', 'expDuration','meanLUM','LMean','targSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','randCond','rgbRatios_raw','rgbRatios'); % save parameters

catch exception
    if ~isempty(dev)
        TouchQueueRelease(dev);
    end
    Screen('CloseAll');
    disp(exception.message)
    for i = 1:numel(exception.stack)
        exception.stack(i)
    end

    dataFile=fullfile( sprintf('%s_AIM_CS_IncremFailedRun_.mat', testSName));  % file path to unique Mat files

    save(dataFile,'trialRecord','expDuration','meanLUM','targSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect', 'examtime'); % save parameters
end
% end

% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Updater function
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject ] = fitPFuncnAFC_AIM_ColorDetect( level, response, nAFC, graphData )
% Fit data from nAFC using cumulative gaussian psychometric function.
% Function returns fitted parameters, goodness of fit, and 95%
% confidence intervals for the fitted parameters.
% % example data
% level=[0.500000000000000;0.445625469066873;0.397164117362141;0.353972892192069;0.315478672240097;0.281170662595175;0.250593616813636;0.223341796075482;0.199053585276749;0.177406694616788;0.158113883008419;0.140919146563223;0.125594321575479;0.111936056928417;0.0997631157484441;0.125594321575479;0.111936056928417;0.0997631157484441;0.0889139705019463;0.0792446596230558;0.0997631157484441;0.0889139705019463;0.0792446596230558;0.0706268772311378;0.0629462705897085;0.0561009227150983;0.0500000000000001;0.0445625469066874;0.0397164117362142;0.0353972892192070;0.0445625469066874;0.0397164117362142;0.0353972892192070;0.0315478672240097;0.0397164117362142;0.0500000000000001;0.0629462705897085;0.0561009227150983;0.0706268772311379;0.0629462705897085;0.0792446596230559;0.0997631157484443;0.0889139705019464;0.0792446596230559;0.0706268772311380;0.0629462705897086;0.0561009227150984;0.0706268772311380;0.0889139705019465;0.111936056928417];
% response=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;0;1;1;1;1;1;1;1;1;1;0;1;1;1;0;0;0;1;0;1;0;0;1;1;1;1;1;0;0;0;0];
% nAFC=4;
% graphData=1;
if nAFC>0; pGuess=1/nAFC; % guess rate
else pGuess=0;
end
h=figure('visible','off');
arrangedData=arrangeData(single(level), response);
tLevel=arrangedData(:,1);
nTests=arrangedData(:,2);
pCorrect=arrangedData(:,4);
dataErrEst=arrangedData(:,8);
%   nAFC cumulative gaussian
ft = fittype( @(thresh, slope,x)(pGuess + (1-pGuess).*(0.5+0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst);
% [fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst,  'Upper',[max(tLevel) 0.21],  'Lower',[min(tLevel) min(diff(tLevel))]);
%    [fitobject, gof] = fit(tLevel, pCorrect,ft);%, 'StartPoint', [mean(tLevel) std(tLevel)]);
thresholdEst=fitobject.thresh;
ci = confint(fitobject);
expVals=fitobject(tLevel); % evaluate fit at test levels
chiSqFitReject=chi2gof(sum(((expVals-pCorrect).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tLevel), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)
if graphData==1
    if any(isnan(ci))
        ci95=zeros(size(tLevel));
    else
        ci95 = predint(fitobject,tLevel,0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    end
    %              scatter(tLevel, pCorrect,nTests+20, 'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    %              hold on
    errorbar(tLevel, pCorrect, dataErrEst, dataErrEst, 'bo');%, 'MarkerSize',nTests+2);
    hold on
    ylim([0,1]);
    %     xlim([min(tLevel)-0.75,max(tLevel)+0.75]);
    plot(fitobject, 'r-');
    if ~any(isnan(ci95(:))) % no NaNs in confidence intervals
        plot(tLevel,ci95,'g--'); % plot 95% confidence interval on fit line
    end
    hold off
    legend HIDE
    box off
    %     set(gca, 'XScale', 'log');
    xlabel('Test Level');
    ylabel('Proportion Correct');
    title(sprintf('Threshold=%.4f, 95%% CI (%.4f,%.4f)', thresholdEst, ci(1,1), ci(2,1)));
end
end


%% AIM Color detect error function
%Version 02/24/2023,Jsk
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject, AIMcompared] = fitPFuncOriErr_AIMColDetect( tOri, mErr, ranErr, graphData,eye,sName, criterion  )
% Fit data from angular error (target - indicated responses) using cumulative gaussian psychometric function.
%% Input
% tOri= tested orientation angles
% mErr= measured error angle
% ranErr= max error angle ( C=90 in each direction and gratings/bipartie  section =45 )
% graphData= plot figure with outcomes
% eye= which eye has been tested? OS, OD, OU
% sName= subject name
% criterion= which point at the function should be choosen? depending on criterion e.g. 60% correct== 72 degree error in each direction
% may vary between AIM tests, AIMtype will indicated which fixed variables should be choosen
%% Output
% fitobject= Function returns fitted parameters for semi and constraint fits
% gof= goodness of fit for semi and constraint fits
% thresholdEst=threshold estimates for semi and constraint fits
% ci=95 % Confidence intervals for the fitted parameters for semi and constraint fits
% chiSqFitReject = to test whether data are indeed signifcantly different from zero for semi and constraint fits
% AIMcompared= Contains AIM semi constraint and its comparison e.g. for visual acuity AIM VA and ETDRS equivalent (fully constraint)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% default input values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    criterion=60; % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) AIM parameters: fit semi-constraint function  i.e. boundaries constraint
ft = fittype( @(thresh, slope, minErr, x)(minErr + (ranErr-minErr).*(0.5-0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
UBslopeconstraint=0.15; % keep it from getting too shallow
LBslopeconstraint=0.05; % keep it from getting too steep
[fitobject.semiconstraint, gof.semiconstraint] = fit(tOri, mErr,ft, 'StartPoint', [mean(tOri) std(tOri) min(mErr)], 'Upper',[max(tOri) UBslopeconstraint ranErr*0.5],  'Lower',[min(tOri) LBslopeconstraint 1]); % semi-constraint
thresholdEst.semiconstraint=fitobject.semiconstraint.thresh;
ci.semiconstraint = confint(fitobject.semiconstraint);
% Find X value were minErr
criterion_minErr= 0.01;% 1 percent of max error +min Err
AIMcompared.semiconstraint.X_minErrExceeded= (fitobject.semiconstraint.thresh)+(sqrt(2)*fitobject.semiconstraint.slope)*(erfinv(criterion_minErr-0.5)/(-0.5));
AIMcompared.semiconstraint.ROSDI=range([thresholdEst.semiconstraint AIMcompared.semiconstraint.X_minErrExceeded]); % calculate the range from noise to threshold (VA) value in logMAR
%compare to semi-constraint ETDRS equivalent
expVals=fitobject.semiconstraint(tOri); % evaluate fit at test levels
chiSqFitReject.semiconstraint=chi2gof(sum(((expVals-mErr).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tOri), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)
AIMcompared.semiconstraint.criterionpoint=(sqrt(2).*fitobject.semiconstraint.slope)*(erfinv(((criterion-fitobject.semiconstraint.minErr )/((ranErr-fitobject.semiconstraint.minErr))-0.5)/-0.5))+thresholdEst.semiconstraint; % ang error func

%% 2) fit unconstraint function with 1 free parameters, i.e. threshold while slope and min error are fixed
ft = fittype( @(thresh, slope, minErr, x)(minErr + (ranErr-minErr).*(0.5-0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
[fitobject.unconstraint, gof.unconstraint]= fit(tOri, mErr,ft,'StartPoint', [mean(tOri) std(tOri) min(mErr)]);% unconstraint
thresholdEst.unconstraint=fitobject.unconstraint.thresh;
ci.unconstraint = confint(fitobject.unconstraint);
% Find X value were minErr
criterion_minErr= 0.01;% 1 percent of max error +min Err
AIMcompared.unconstraint.X_minErrExceeded= (fitobject.unconstraint.thresh)+(sqrt(2)*fitobject.unconstraint.slope)*(erfinv(criterion_minErr-0.5)/(-0.5));
AIMcompared.unconstraint.ROSDI=range([thresholdEst.unconstraint AIMcompared.unconstraint.X_minErrExceeded]); % calculate the range from noise to threshold (VA) value in logMAR
expVals=fitobject.unconstraint(tOri); % evaluate fit at test levels
chiSqFitReject.unconstraint=chi2gof(sum(((expVals-mErr).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tOri), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)
AIMcompared.unconstraint.criterionpoint=(sqrt(2).*fitobject.unconstraint.slope)*(erfinv(((criterion-fitobject.unconstraint.minErr )/((ranErr-fitobject.unconstraint.minErr))-0.5)/-0.5))+thresholdEst.unconstraint; % ang error func

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphData==1
    %figure()
    %% plot AIM semiconstraint fit
    tOriSort=sort(tOri);
    ci95 = predint(fitobject.unconstraint,tOriSort, 0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    scatter(tOri, mErr);
    hold on
    plot(fitobject.unconstraint, '-');
    h.LineWidth=1;
    h.Color=[0.9290, 0.6940, 0.1250];
    if ~any(isnan(ci95(:))) % no NaNs in confidence intervals
        plot(tOriSort,ci95','--', 'Color',[0, 0.4470, 0.7410],'LineWidth',0.25); %  95% CIs
    end
    hold off
    legend HIDE
    box off
    xlabel('Linear Contrast'); % label the x axis;JH11/18/2022:was log10(Contrast), otherwise can not use error function
    ylabel('Angular Error [ \circ ]'); % label y axis
    title(sprintf('%s  | %.2fdeg |threshold %.2f%% |CI(%.2f %.2f)', sName, 100*10.^thresholdEst.unconstraint, 100*10.^ci95(1), 100*10.^ci95(2)));
    ylim([0 180]); % fix upper and lower bounds

    %subtitle([' Noise= ' num2str(fitobject.semiconstraint.minErr,2), '\circ ','| ROSDI= ',num2str(AIMcompared.semiconstraint.ROSDI,2)]);
end
end

