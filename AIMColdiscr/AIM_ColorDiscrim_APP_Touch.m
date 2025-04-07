% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/20/2023

function [output_graph, testSName, trialRecord, randCond, expDuration, correction, angleFitobject, AIMcompared, confidenceInterval95] = AIM_ColorDiscrim_APP_Touch(view_distance, Patient_Name, eye, number_of_trials, correction_type)

% clear all; close all; commandwindow; % close windows and type into command window
% Patient_Name = 'nj'; testEye='3'; Patient_ID='bhjd124'; view_distance= '40';correction_type = 1; number_of_trials= 2; eye=1;

fprintf("#######################AIM_Color_Discrimination Starting_%s#######################\n", datestr(now,'mm-dd-yyyy HH-MM-SS'))
examtime = datetime;
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
% prompt = {'Subject Initials', 'Pedestal Contrast', 'C Size (deg)','Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)','With(1) of Without(0) Correction','Which eye? left=1,right=2,both=3'};
% dlg_title = 'AIM Color Discrim';
% num_lines = 1;

sName=char(Patient_Name); % assign experimenter values to experiment
pedContrast=str2num(char('0.2')); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
targSizeDeg=str2num(char('2')); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
cellSizeDeg=str2num(char('5')); % diameter of cell rings
nRows=str2num(char('4')); % $ rows and columns on each chart
nCols=str2num(char('4'));
nTrials=number_of_trials; % # charts to run

[width, height]=Screen('DisplaySize', whichScreen); % adaptive screen width measurement in mm
scrnWidthCm = width/10; % converting mm to cm

viewDistCm=str2num(char(view_distance)); % observer's distance from screen
correction=correction_type; %with (1) or without (2) spectacle correction
testEye=eye; %1=left,2=right, 3=both

meanLUM=0.5;
rng('default'); % seed random number generator
gammaVal=2.22; % gamma on this system
LMean=round(255*meanLUM.^(1/gammaVal)); % Gamma-corrected background
noiseRateHz=14; % update rate of noise pedestal
ringLineDeg=str2num(char('0.2')); % indication ring width in degree
nAFC=8; % angle for error scoring
saveImage=0; % save stimuli or not
interCellGapDeg=1; %changed from
ranErrorDeg=45;
minCRange=0.002; % minium range of different angle %jh 11/16/2022: change min contrast range from 0.01 to 0.004

rng('shuffle')
% load('equilumColorTable_thresh001step2deg');
% equilumColorTable=equilumColorTable_thresh;

load('equilumColorTable_new.mat');
stepsize=0.001; % color lookup table stepsize
threshMultiple = 6;
nColors=4;% nColors=length(hsvAngs); %size(rgbRatios,1);% how many colors to test %jh: needs to replace color space

colAngs_raw=linspace(0,2*pi,8+1); % get angle of n axes in rad,has not rounded to resolution limit

colAngs_raw = colAngs_raw(1:2:end);
%purple S+
% randCond = randperm(nColors);

randCond = randperm(4);

colAngs=colAngs_raw(randCond);% randomization of rows of color directions % order of colors: purple, red, yellow, green
% randCond = [1,2,3,4,5,6,7,8];colAngs=colAngs_raw;
ndigit = 3; % round to nth digit based on color lookup table step size, 0.001step size round to 3th digit

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

dataFile=fullfile(cd,'Patient Data', sprintf('%s_AIM_ColDiscrim.mat', testSName));  % file path to unique Mat fes

dev = min(GetTouchDeviceIndices([], 1));
if ~isempty(dev)
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
end

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    screenSize = get(0,'screensize');
    Screen('Preference','VisualDebugLevel', 0);  % Remove pysch toolbox warning (Red screen
    % [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean,screenSize);
    [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean); %,[0 0 3839, 2160]   ,[0 0 2800 1800]);


    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels

    global pixPerDeg
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree

    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    responseRingPenWidth=pixPerDeg*ringLineDeg; % line width for response ring
    responseRingRect=[0 0 cellDiamPix-responseRingPenWidth/2 cellDiamPix-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);
    testGapWidthDeg=360/(5*pi); % gap is always 22.9 deg (1/5 of the line width)

    targSizePix=round(targSizeDeg*pixPerDeg); % size of target on this system
    stimPix=targSizePix; % targSizePix or cellDiamPix
    rndWin=makeDisk(cellDiamPix);
    stimWin=makeDisk(cellDiamPix,stimPix/2);

    %     myEdge=ones(stimPix);
    %     myEdge(:, 1:stimPix/2)=0; % left half 0
    %     myEdge=imgaussfilt(myEdge,pixPerDeg*0.15);
    mystimPix=makeStim(cellDiamPix, stimPix);% make target of required size in noise background size
    colLayer=zeros(cellDiamPix);

    checkSizePix=4;% always 0.25 deg: round(0.25*pixPerDeg); % size of random noise pedestal checks
    nPrecomputedFrames=14;%14 noise pages and 14 samples / sec; JSk, JH, 01/30/2023
    srcRect=[0 0 cellDiamPix cellDiamPix]; % rect for this sized cell

    nextSizeDeg=4; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left

    responseRingCol=[LMean LMean LMean]; % JH,Jsk- changed response ring to mean lum
    responseWidthDeg=10;

    if ~isempty(dev)
        TouchQueueCreate(windowPtr, dev);
        TouchQueueStart(dev);
        KbReleaseWait;
        HideCursor(windowPtr);
    else
        ShowCursor('Arrow'); % show arrow response cursor
    end

    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,100);  % put some instructions on screen
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
    else
        [mx,my,buttons] = GetMouse; % wait for mouse button release before starting, experiment
        while any(buttons) % if already down, wait for release
            [mx,my,buttons] = GetMouse(windowPtr);
        end
        while ~any(buttons) % wait for new press
            [mx,my,buttons] = GetMouse();
        end
    end

    instructions = 'Each cell contains two colors.\n Tap where you see the boundary for each cell. \n If you cannot see the boundary, make a guess';
    instructImage = imread("AIM_ColorDiscrim_Image_rotate.jpg");
    Texture = Screen('MakeTexture', windowPtr, instructImage);
    rect = [1050, 1200, 1900, 1600];

    Screen('DrawTexture', windowPtr, Texture, [], rect);

    [~, ~, textBox] = DrawFormattedText(windowPtr,  instructions, 'center', 'center',0, LMean); % give instructions to subject
    Screen('Flip', windowPtr); % flip to the information screen
    WaitSecs(1); % UX consideration: provide time to prepare to start
    command = audioread('female_AIMColDiscrim_instruction.mp3'); % sound file when mouse is clicked
    sound(command, 22000);
    WaitSecs(8);
    %     Screen('CloseAll'); % close all windows

    Screen('TextSize',windowPtr,48);


    % slopeEsts=2*ones(size(targSizeDeg)); % start with a rough estimate of slope
    % minErrEsts=2*ones(size(targSizeDeg)); % start with a rough estimate of internal error
    % startContrastBnds = [2*stepsize, 10*stepsize, 20*stepsize, 10*stepsize, 2*stepsize, 10*stepsize, 20*stepsize, 20*stepsize;...
    %     0.8*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors),0.8*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors)]';% starting index in deg or rad for different colors
    % startContrastBnds = [2*stepsize, 10*stepsize, 20*stepsize, 10*stepsize;...
    % 0.8*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors)]';% starting index in deg or rad for different colors

    %     startContrastBnds = [2*stepsize, 10*stepsize, 20*stepsize, 10*stepsize;...
    %         0.8*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors),2*(2*pi/nColors)]';% starting index in deg or rad for different colors
    startContrastBnds = [2*stepsize, 10*stepsize, 2*stepsize, 10*stepsize; % chnging min sizes from above to stesize levels
        0.5*(pi/nColors), 1.25*(2*pi/nColors), 0.5*(pi/nColors), 1.25*(2*pi/nColors)]'; % iterating between angles

    startContrastBnds=startContrastBnds(randCond,:); % Randomizing color boundary order
    startContrastBnds = round(startContrastBnds,3); % rounding color boundary to 3 digits
    %      for ii=1:length(randCond);startContrastBnds(ii,:)=startContrastBnds(randCond(ii),:);end % re-order starting index to match color order
    % the max angle around each test axis is the angle between the 2
    % adjacent angles, ie 2*(2pi/nColors). the min angle is the resolution
    % of the angle~RGB lookup table, which is 0.01 in rad.

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
    chart =0;

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
                angErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 45 deg
                %                 fitObj=fitPFuncOriErr(testLevel',abs(angErr)',1); % fit with current orientation errors
                fitObj=fitPFuncnAFC_AIM_ColorDiscrim(testLevel, respCorrect, nAFC, 0); % use matlab's fit function to fit cumulative Gaussian to the data
                csBnds = confint(fitObj, 0.95); % 95% confidence intervals on contrast threshold estimate
                if any(isnan(csBnds)) % fit failed to estimate acuity and 95% CIs
                    contrastBound(1)=fitObj.thresh-2*fitObj.slope; % lower bound is threshold - 2 std
                    contrastBound(2)=fitObj.thresh+2*fitObj.slope; % upper bound is threshold + 2 std
                else
                    contrastBound(1)=csBnds(1);
                    contrastBound(2)=csBnds(2);
                end
                if any(isnan(contrastBound)) % fit failed to estimate slope
                    contrastBound=startContrastBnds(colNum,:); % use defaults
                end
                % jh: default [min max] should be [2*stepsize 2*(2*pi/nColors)];
                minTestLevel=max([2*stepsize contrastBound(1)]); % jh value adjusted
                minTestLevel=min([2*(2*pi/nColors)-minCRange minTestLevel]); % prevent minTestLevel from getting too large (>0)
                % when fitted curve is flat,abs(contrastBound) becomes too large, so force it to be log10(1-pedContrast)-minCRange
                maxTestLevel=min([2*(2*pi/nColors) contrastBound(2)]); % uppper 95% CI in arcmin and maximum size within target circle
            else
                contrastBound=startContrastBnds(colNum,:); % first trial: use defaults
                minTestLevel=contrastBound(1);
                maxTestLevel=contrastBound(2);
            end

            if maxTestLevel<minTestLevel
                maxTestLevel=2*(2*pi/nColors);
                minTestLevel=2*(2*pi/nColors)-minCRange;
            end % prevent maxTestLevel from going outside of the range
            cRange=maxTestLevel-minTestLevel;
            % if cRange<minCRange % if range is too small, then increase
            %     midLevel=minTestLevel+(maxTestLevel-minTestLevel)/2;
            %     minTestLevel=max([startContrastBnds(colNum,1) midLevel-minCRange/2]);
            %     maxTestLevel=min([startContrastBnds(colNum,2) midLevel+minCRange/2]);
            % end
            vbl= Screen('Flip', windowPtr); % clear screen
            buttonPressTime=vbl;
            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(trialNo, colNum).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(trialNo, colNum).trialSeed); % seed the random number generator

            trialRecord(trialNo,colNum).targLocs=Shuffle(1:nCells); % pick random target locations for target contrasts
            trialRecord(trialNo,colNum).targContrast=linspace(minTestLevel, maxTestLevel, nCells); % linear spaced range of color angle between upper and lower test bounds
            trialRecord(trialNo,colNum).targContrastPerLoc(trialRecord(trialNo,colNum).targLocs)=trialRecord(trialNo,colNum).targContrast; % fill target locations with corresponding contrast

            trialRecord(trialNo,colNum).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,colNum).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,colNum).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,colNum).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

            for cellNum=1:nCells % draw all random orientation rings
                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord(trialNo,colNum).targOri(cellNum)=testAngle; % update record for test angle

                % JH modified color space
                color1=colAngs(colNum)-trialRecord(trialNo,colNum).targContrastPerLoc(cellNum)/2;
                if color1<0;
                    color1=round(color1+round(2*pi,ndigit),ndigit);
                elseif color1>round(2*pi,ndigit);
                    color1=round(color1-round(2*pi,ndigit),ndigit);
                else color1=round(color1,ndigit);
                end
                color2=colAngs(colNum)+trialRecord(trialNo,colNum).targContrastPerLoc(cellNum)/2;
                if color2<0; color2=round(color2+round(2*pi,ndigit),ndigit);elseif color2>round(2*pi,ndigit); color2=round(color2-round(2*pi,ndigit),ndigit);else color2=round(color2,ndigit);end
                RGBvals1=equilumColorTable(equilumColorTable(:,1)==color1,2:4);% RGB in [0 1]
                RGBvals2=equilumColorTable(equilumColorTable(:,1)==color2,2:4);% RGB in [0 1]
                RGBvals1=threshMultiple*((RGBvals1*2)-1);% RGB+pedastal scale to [-1 1]
                RGBvals2=threshMultiple*((RGBvals2*2)-1);% RGB+pedastal scale to [-1 1]

                myTarg=imrotate(mystimPix, -testAngle-90,'nearest', 'crop'); % rotate the C by required amount
                myTarg= 0.5*myTarg+0.5;

                for preFrameNum=1:nPrecomputedFrames % make as manytextures as requested
                    myNoiseRGB=zeros([cellDiamPix cellDiamPix 3]); % blank stimulus
                    myNoise=randn(round(cellDiamPix/checkSizePix)); % random downsampled noise
                    if checkSizePix>1 % need to resize the noise
                        myNoise=imresize(myNoise, [cellDiamPix cellDiamPix], 'nearest');
                    end
                    myNoise=myNoise/max(abs(myNoise(:))); % scale -1 to +1
                    myNoise=pedContrast*(myNoise.*rndWin); % scale noise to pedestal level--circular noise
                    for rgbLayer=1:3
                        colLayer=myTarg.*RGBvals1(rgbLayer)+(1-myTarg).*RGBvals2(rgbLayer);
                        %colLayer(rndWin==0)=0;
                        colLayer(stimWin==0)=0; % restrict target area
                        myNoiseRGB(:,:,rgbLayer)=myNoise+colLayer;
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

            %[~, ~, textBox] = DrawFormattedText(windowPtr, 'Click on the orientation of each C:','center', 100, 0); % give instructions to subject

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
            frameNum=0;

            while clickedExit==0 % keep checking until subject click exit location
                frameNum=frameNum+1;
                %                 Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
                %                 DrawFormattedText(windowPtr, sprintf('%s', 'next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
                for cellNum=1:nCells % work through each cell
                    Screen('DrawTexture', windowPtr, myTex(cellNum,1+mod(frameNum,nPrecomputedFrames)), srcRect,destRect(cellNum,:)); % draw texture to on screen location
                    Screen('FrameOval',windowPtr,responseRingCol*0.8, destRect(cellNum,:), responseRingPenWidth);  % draw box around cells
                    if trialRecord(trialNo,colNum).stimSeen(cellNum)==1 % this cell clicked - put a ring on it
                        Screen('FrameOval',windowPtr,[LMean*0.9 LMean*0.9 LMean*0.9], destRect(cellNum,:), responseRingPenWidth);  % draw ring around cells
                        Screen('FrameArc',windowPtr,[LMean*1.2 LMean*1.2  LMean*1.2 ],destRect(cellNum,:), 90+trialRecord(trialNo,colNum).matchOri(cellNum)-responseWidthDeg/2,responseWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc - add 90 because PTB 0 = up
                        Screen('FrameArc',windowPtr,[LMean*1.2  LMean*1.2  LMean*1.2 ],destRect(cellNum,:),-90+trialRecord(trialNo,colNum).matchOri(cellNum)-responseWidthDeg/2,responseWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc - add 90 because PTB 0 = up
                    end
                end
                vbl= Screen('Flip', windowPtr, vbl+(1/noiseRateHz)); %i.e. 1/8 Hz==125ms; draws noise frame every 125ms

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
                            if evt.MappedX < destRect(16, 3) && evt.MappedX > destRect(1,1) && evt.MappedY < destRect(16, 4) && evt.MappedY > destRect(1,2)
                                beepData = audioread('pop.wav'); % sound file when mouse is clicked
                                sound(beepData, 50000);
                            end
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

                                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                                    [~,respNum]=min(mouseDistFromEachBox(:));
                                    %                                     if any(buttons) % wait 250 msec between clicks
                                    buttonPressTime=vbl; % note screen time when button was pressed
                                    if mx > nextRect(1) && mx < nextRect(3) && my >nextRect(2) && my< nextRect(4) && sum(trialRecord(trialNo,colNum).stimSeen)==nCells % observer clicked finished word
                                        clickedExit=1;
                                        beepData = audioread('good.wav'); % sound file when mouse is clicked
                                        sound(beepData, 50000);
                                    else
                                        trialRecord(trialNo,colNum).matchOri(respNum)=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                                        respErr(1)=diff(unwrap([trialRecord(trialNo,colNum).targOri(respNum),trialRecord(trialNo,colNum).matchOri(respNum)]/180*pi)*180/pi); % error of 0 and 180 deg because grating flips at 180
                                        respErr(2)=diff(unwrap([trialRecord(trialNo,colNum).targOri(respNum),trialRecord(trialNo,colNum).matchOri(respNum)+180]/180*pi)*180/pi);
                                        [~,minErrLoc]=min(abs(respErr)); % observer's error is closer of 2
                                        trialRecord(trialNo,colNum).oriErr(respNum)=respErr(minErrLoc); % calculate difference betwen actual and reported orientation;% calculate difference betwen actual and reporte orientation
                                        trialRecord(trialNo,colNum).stimSeen(respNum)=1;
                                        if abs(trialRecord(trialNo,colNum).oriErr(respNum))<(180/nAFC) % if within nAFC error
                                            trialRecord(trialNo,colNum).respCorrect(respNum)=1; % score correct
                                        else
                                            trialRecord(trialNo,colNum).respCorrect(respNum)=0; % score incorrect
                                        end

                                    end
                                    %  end
                                    trialRecord(trialNo,colNum).chartTime=toc(chartStart);
                                else
                                    blobcol{id} = [];
                                end
                            end
                        end
                    end
                    if sum( trialRecord(trialNo,colNum).stimSeen)==nCells && nTrials==trialNo && colNum==nColors%change color of next button to green once all cells have been clicked, say NEXT
                        Screen('FillOval',windowPtr,[LMean*0.9 LMean*0.9 LMean*0.9], nextRect);  % draw exit button, green
                        DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [], [],[],[],[],[],nextRect); %  char(26)
                    elseif sum( trialRecord(trialNo,colNum).stimSeen)==nCells
                        Screen('FillOval',windowPtr,[LMean*0.9 LMean*0.9 LMean*0.9], nextRect);  % draw exit button, green
                        DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [], [],[],[],[],[],nextRect); %  char(26)
                    end
                end
            end
            trialRecord(trialNo,colNum).chartTime=toc(chartStart);
            if saveImage==1 % saving demo image?
                myImfileName=sprintf('AIM_ColDiscrim%d.jpg', colNum); % create filename
                myImage=Screen('GetImage', windowPtr); % grab screnshot
                imwrite(myImage,myImfileName); % write screenshot to image file
            end
            %% Were the responses significantly different from random chances? If so continue, otherwise show the following screen
            temp_y_90(1:length(trialRecord(trialNo,colNum).oriErr(:)))=ranErrorDeg;
            [h,p] =  ttest(abs(trialRecord(trialNo,colNum).oriErr(:)),temp_y_90');

            if h==0 %  minimum response error is not significantly different from random error either 45 or 90 deg
                % Note: random responses, not different from 90/45
                if chart < (nTrials * nColors)
                    instcommand1 = audioread('female_AIMColDiscrim_warning.mp3'); % sound file when mouse is clicked
                    sound(instcommand1, 22000)
                    WaitSecs(7)
                end
                trialRecord(trialNo,colNum).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response

            elseif h==1 %minimum response error is significantly different from random error either 45 or 90 deg
                % min error is significantly different but we are not sure
                % if its larger than 45/90 deg therefore we do t test
                [h1,p1] =  ttest(abs(trialRecord(trialNo,colNum).oriErr(:)),ranErrorDeg, 'Tail','right');
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
                        instcommand1 = audioread('female_AIMColDiscrim_warning.mp3'); % sound file when mouse is clicked
                        sound(instcommand1, 22000)
                        WaitSecs(7)
                    end
                    trialRecord(trialNo,colNum).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response
                end
            end
        end % end SFs loop
    end % end trials loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    if ~isempty(dev)
        TouchQueueStop(dev);
        TouchQueueRelease(dev);
        ShowCursor(windowPtr);
    end
    Screen('CloseAll'); % close all windows
    save(dataFile,'examtime','correction','testEye','trialRecord', 'expDuration','meanLUM','targSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','threshMultiple','equilumColorTable','randCond','colAngs_raw','colAngs'); % save parameters
    output_graph(1) = figure('visible','off');
    [~,idx]=sort(randCond);
    for colNum=1:nColors % plot figure with fits for interleaved conditions - attempt to fit all data
        %subplot(ceil(nColors^0.5),ceil(nColors/ceil(nColors^0.5)),colNum);
        % if idx(colNum)==1 %plot colAngs_raw value 1==0
        subplot(1,nColors,colNum); % create new figure
        testLevel=[]; % set blank list of all levels for this SF
        respErr=[]; % blank list of all stimuli seen for this SF
        stimSeen=[]; % blank list of all stimuli seen for this SF
        for trialSoFar=1:nTrials % work through all trials
            testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
            respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
            stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
        end
        respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
        %             testLevel=10.^testLevel;% JH added 11/18/2022
        [angleFitobject{colNum},~,angThreshCS{colNum},~,~,AIMcompared{colNum}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data

        if colNum==1
            confidenceInterval95.purple = confint(angleFitobject{colNum}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
            fitLB(colNum)=confidenceInterval95.purple(1,1);
            fitUB(colNum)=confidenceInterval95.purple(2,1);
            subtitle('Purple|S+');
        elseif colNum==2
            confidenceInterval95.red = confint(angleFitobject{colNum}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
            fitLB(colNum)=confidenceInterval95.red(1,1);
            fitUB(colNum)=confidenceInterval95.red(2,1);
            subtitle('Red|L-M');
        elseif colNum==3
            confidenceInterval95.yellow = confint(angleFitobject{colNum}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
            fitLB(colNum)=confidenceInterval95.yellow(1,1);
            fitUB(colNum)=confidenceInterval95.yellow(2,1);
            subtitle('Yellow|S-')
        elseif colNum==4
            confidenceInterval95.green = confint(angleFitobject{colNum}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
            fitLB(colNum)=confidenceInterval95.green(1,1);
            fitUB(colNum)=confidenceInterval95.green(2,1);
            subtitle('Green|M-L')
        end
    end
    %         if idx(colNum)==1%plot colAngs_raw value 1==
    %             subplot(2,nColors/2,1);
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
    %             %             testLevel=10.^testLevel;% JH added 11/18/2022
    %             [angleFitobject{1},~,angThreshCS{1},~,~,AIMcompared{1}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    %             confidenceInterval95= confint(angleFitobject{idx(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    %             fitLB(1)=confidenceInterval95(1,2);
    %             fitUB(1)=confidenceInterval95(2,2);
    %             subtitle(['Purple-Red'])

    %         elseif idx(colNum)==2%plot colAngs_raw value 1==0
    %             subplot(2,nColors/2,2);
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
    %             %             testLevel=10.^testLevel;% JH added 11/18/2022
    %             [angleFitobject{2},~,angThreshCS{2},~,~,AIMcompared{2}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    %             confidenceInterval95.red = confint(angleFitobject{idx(colNum)}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    %             fitLB(2)=confidenceInterval95.red(1,2);
    %             fitUB(2)=confidenceInterval95.red(2,2);
    %
    %             subtitle(['Red|L-M'])
    %
    % %         elseif idx(colNum)==2%plot colAngs_raw value 1==0
    % %             subplot(2,nColors/2,2);
    % %             testLevel=[]; % set blank list of all levels for this SF
    % %             respErr=[]; % blank list of all stimuli seen for this SF
    % %             stimSeen=[]; % blank list of all stimuli seen for this SF
    % %             for trialSoFar=1:nTrials % work through all trials
    % %                 testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
    % %                 respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
    % %                 stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
    % %             end
    % %             respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
    % %             %             testLevel=10.^testLevel;% JH added 11/18/2022
    % %             [angleFitobject{2},~,angThreshCS{2},~,~,AIMcompared{2}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    % %             confidenceInterval95= confint(angleFitobject{idx(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    % %             fitLB(2)=confidenceInterval95(1,2);
    % %             fitUB(2)=confidenceInterval95(2,2);
    % %             subtitle(['Red-Yellow'])
    %
    %         elseif idx(colNum)==3%plot colAngs_raw value 1==0
    %             subplot(2,nColors/2,3);
    %             testLevel=[]; % set blank list of all levels for this SF
    %             respErr=[]; % blank list of all stimuli seen for this SF
    %             stimSeen=[]; % blank list of all stimuli seen for this SF
    %             for trialSoFar=1:nTrials % work through all trials
    %                 testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
    %                 respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
    %                 stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
    %             end
    %             respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
    %             %testLevel=10.^testLevel;% JH added 11/18/2022
    %             [angleFitobject{3},~,angThreshCS{3},~,~,AIMcompared{3}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    %             confidenceInterval95.yellow = confint(angleFitobject{idx(colNum)}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    %             fitLB(3)=confidenceInterval95.yellow(1,2);
    %             fitUB(3)=confidenceInterval95.yellow(2,2);
    %             subtitle(['Yellow|S-'])
    %
    % %         elseif idx(colNum)==3%plot colAngs_raw value 1==0
    % %             subplot(2,nColors/2,3);
    % %
    % %             testLevel=[]; % set blank list of all levels for this SF
    % %             respErr=[]; % blank list of all stimuli seen for this SF
    % %             stimSeen=[]; % blank list of all stimuli seen for this SF
    % %             for trialSoFar=1:nTrials % work through all trials
    % %                 testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
    % %                 respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
    % %                 stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
    % %             end
    % %             respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
    % %             % testLevel=10.^testLevel;% JH added 11/18/2022
    % %             [angleFitobject{3},~,angThreshCS{3},~,~,AIMcompared{3}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    % %             confidenceInterval95= confint(angleFitobject{idx(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    % %             fitLB(3)=confidenceInterval95(1,2);
    % %             fitUB(3)=confidenceInterval95(2,2);
    % %             subtitle(['Yellow-Green'])
    %
    %         elseif idx(colNum)==4%plot colAngs_raw value 1==0
    %             subplot(2,nColors/2,4);
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
    %             %             testLevel=10.^testLevel;% JH added 11/18/2022
    %             [angleFitobject{4},~,angThreshCS{4},~,~,AIMcompared{4}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName);
    %             confidenceInterval95.green = confint(angleFitobject{idx(colNum)}.semiconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    %             fitLB(4)=confidenceInterval95.green(1,2);
    %             fitUB(4)=confidenceInterval95.green(2,2);
    %             subtitle(['Green|M-L'])
    %
    %
    % %         elseif idx(colNum)==4%plot colAngs_raw value 1==0
    % %             subplot(2,nColors/2,4);
    % %
    % %             testLevel=[]; % set blank list of all levels for this SF
    % %             respErr=[]; % blank list of all stimuli seen for this SF
    % %             stimSeen=[]; % blank list of all stimuli seen for this SF
    % %             for trialSoFar=1:nTrials % work through all trials
    % %                 testLevel=[testLevel trialRecord(trialSoFar,idx(colNum)).targContrastPerLoc]; % concatenate all test sizes
    % %                 respErr=[respErr trialRecord(trialSoFar,idx(colNum)).oriErr]; % concatenate all orientation errors
    % %                 stimSeen=[stimSeen trialRecord(trialSoFar,idx(colNum)).stimSeen]; % concatenate all responses made
    % %             end
    % %             respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
    % %             %testLevel=10.^testLevel;% JH added 11/18/2022
    % %             [angleFitobject{4},~,angThreshCS{4},~,~,AIMcompared{4}]=fitPFuncOriErr_AIMColDiscrim(rad2deg(testLevel'),abs(respErr)',ranErrorDeg,1,testEye,sName); % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
    % %             confidenceInterval95= confint(angleFitobject{idx(colNum)}.unconstraint, 0.95); % 95% confidence intervals on contrast threshold estimate
    % %             fitLB(4)=confidenceInterval95(1,2);
    % %             fitUB(4)=confidenceInterval95(2,2);
    % %             subtitle(['Green-Purple'])

    %         end
    %     end
    save(dataFile,'examtime','correction','testEye','colAngs','trialRecord','angleFitobject','AIMcompared', 'angThreshCS', 'expDuration','meanLUM','targSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','threshMultiple','equilumColorTable','randCond','colAngs_raw','colAngs', 'confidenceInterval95','fitLB', 'fitUB'); % save parameters
    %% summary threshold plot
    output_graph(2)=figure('visible','off'); % create new figure;
    for colNum=1:nColors % process data foreach Spatial Frequency
        if colNum==1 %
            e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            e.Marker = 'o';
            e.MarkerSize = 7;
            e.Color = 'k';
            e.CapSize = 15;
        elseif colNum==2 %
            e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            e.Marker = 'o';
            e.MarkerSize = 7;
            e.Color = 'k';
            e.CapSize = 15;
        elseif colNum==3 %
            e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            e.Marker = 'o';
            e.MarkerSize = 7;
            e.Color = 'k';
            e.CapSize = 15;
        elseif (colNum)==4 %
            e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            e.Marker = 'o';
            e.MarkerSize = 7;
            e.Color = 'k';
            e.CapSize = 15;
            %         elseif colNum==5 %
            %             e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            %             e.Marker = 'o';
            %             e.MarkerSize = 7;
            %             e.Color = 'k';
            %             e.CapSize = 15;
            %         elseif colNum==6 %
            %             e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            %             e.Marker = 'o';
            %             e.MarkerSize = 7;
            %             e.Color = 'k';
            %             e.CapSize = 15;
            %         elseif colNum==7 %
            %             e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            %             e.Marker = 'o';
            %             e.MarkerSize = 7;
            %             e.Color = 'k';
            %             e.CapSize = 15;
            %         elseif colNum==8 %
            %             e= errorbar(colNum,angThreshCS{colNum}.semiconstraint, angThreshCS{colNum}.semiconstraint-fitLB(colNum),fitUB(colNum)-angThreshCS{colNum}.semiconstraint);hold on
            %             e.Marker = 'o';
            %             e.MarkerSize = 7;
            %             e.Color = 'k';
            %             e.CapSize = 15;
        end
        title('AIM Color Discrimination','FontSize',20)
        subtitle('Purple|S+    Red|L-M    Yellow|S-    Green|M-L', 'FontSize', 16);
        xlabel('Stimulus types', 'FontSize', 16);
        xlim([0 nColors+1]);
        ylim([0 90]);
        ylabel('Discrimination Angle [\circ]', 'FontSize', 16);
    end
catch exception
    if ~isempty(dev)
        TouchQueueRelease(dev);
    end

    disp(exception.message);
    for i = 1:numel(exception.stack)
        exception.stack(i)
    end
    Screen('CloseAll');

    dataFile=fullfile(cd,'Patient Data', sprintf('%sAIM_ColDiscrimFailedRun.mat', testSName));  % file path to unique Mat fes

    save(dataFile,'trialRecord','expDuration','meanLUM','targSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','threshMultiple','equilumColorTable','randCond','colAngs_raw','colAngs','correction'); % save parameters
end
end


%% make noise and stimulus patch separate
function myStim = makeStim(imSize, targSize)
% Make image of bipartite circular color discrimination stimulus in bacground noise size
% inputs
% imSize: whole size
% targSize: target size
if nargin==1
    targSize=imSize;
end

% myEdge=ones(targSize);
% myEdge(:, 1:targSize/2)=0; % left half 0
% myEdge=imgaussfilt(myEdge,pixPerDeg*0.15);
global pixPerDeg
[X,Y]=meshgrid(linspace(-imSize/2, imSize/2, imSize),linspace(-imSize/2, imSize/2, imSize)); % 2D matrix of radial distances from centre
radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
myStim=zeros(imSize); % set up background

myEdge=ones(targSize);
myEdge(:, 1:targSize/2)=0; % left half 0
myEdge=imgaussfilt(myEdge,pixPerDeg*0.1);% range from 0 to 1
myEdge=myEdge*2-1;% range from -1 to 1

myStim(round((imSize-targSize)/2) +1:round((imSize-targSize)/2) +targSize, round((imSize-targSize)/2) +1:round((imSize-targSize)/2) +targSize)=myEdge; % make disk of 1's of required size
% myLandoltC(abs(X)<imSize/10 & Y<0)=0; % % add vertical top gap
end
% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Updater function
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject ] = fitPFuncnAFC_AIM_ColorDiscrim( level, response, nAFC, graphData )
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
%[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst);
[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst,  'Upper',[max(tLevel) 0.2],  'Lower',[min(tLevel) min(diff(tLevel))]);
%[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst,  'Upper',[max(tLevel) 0.2],  'Lower',[min(tLevel) 0.05]);%Jsk 071323

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

%% AIM Color discrim error function
%Version 02/24/2023,Jsk
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject, AIMcompared] = fitPFuncOriErr_AIMColDiscrim( tOri, mErr, ranErr, graphData,eye,sName, criterion)
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
    criterion=34; % 68degree==62.5%/ 34degree for bi-partie stimulus 4AFC CAD test equivalent criterion Pfit cumulative Gaussian to the ori error data
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) AIM parameters: fit semi-constraint function  i.e. boundaries constraint
ft = fittype( @(thresh, slope, minErr, x)(minErr + (ranErr-minErr).*(0.5-0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
% UBslopeconstraint=0.15; % keep it from getting too shallow - corresponding value for radians
% LBslopeconstraint=0.05; % keep it from getting too steep - corresponding value for radians
UBslopeconstraint=8; % keep it from getting too shallow - corresponding value for degrees
LBslopeconstraint=1; % keep it from getting too steep - corresponding value for degrees
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
    ci95 = predint(fitobject.semiconstraint,tOriSort, 0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    scatter(tOri, mErr);
    hold on
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
    ylabel('Angular Error [\circ]','FontSize',10); % label y axis
    xlim([0,rad2deg(1.6)]);
    ylim([0,ranErr*2]);
    xlabel('Color Angle (radian)','FontSize',10); % label the x axis

    % title_text={sprintf('Threshold: %s', num2str( 360*thresholdEst.unconstraint,3)) ,sprintf('CI: %s %s',num2str(360*ci95(1),3), num2str(360*ci95(2),3))} ;
    % % 'Threshold:' num2str( 360*thresholdEst.unconstraint,3) '|CI:' num2str(360*ci95(1),3) ' ' num2str(360*ci95(2),3)]
    %     title(title_text,'FontSize',10);
    %subtitle([' Noise= ' num2str(fitobject.semiconstraint.minErr,2), '\circ ','| ROSDI= ',num2str(AIMcompared.semiconstraint.ROSDI,2)]);
end
end