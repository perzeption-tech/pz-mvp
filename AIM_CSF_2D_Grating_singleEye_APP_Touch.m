% AIM CSf 2D grating
% Program to measure contrast threshold based on orientation identification of rotated grating C stimulus
% Presents sequence of grids containing random orientation grating stimuli of a range of spatial frequencies and contrasts
% Spatial Frequency range sampled from current acuity estimate to lowest SF possible
% Contrast range of stimuli based on threshold estimate -2 std dev to + 4 stdevs of current slopes
% Observer clicks perceived orientation of each grating and can change response or leave cells empty
% Responses fit with log parabola CSF compbined with cumulative gaussian
% parameters define SF and CS range for next chart
% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/10/2022
function [output_graph, testEye, testSName, csfAcuity, AULCSF, expDuration, fitObj,trialRecord, correction] = AIM_CSF_2D_Grating_singleEye_APP_Touch(view_distance, eye, Patient_Name, Patient_ID, number_of_trials, correction)

% close all; clear all; close all; commandwindow; % close windows and type into command window
% view_distance ='40'; eye=1; Patient_Name='jj'; Patient_ID='sx1'; number_of_trials = 2; correction = 1;

fprintf("#######################AIM_CSF_Starting_%s#######################\n", datestr(now,'mm-dd-yyyy HH-MM-SS'))

examtime= datetime;
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject Initials', 'Background level (0-1)', 'SF range (c/deg)','Parameter Estimate (%)','Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)','With(1) of Without(0) Correction','Which eye? left=1,right=2,both=3','Feedback ON=1, OFF=0'};
dlg_title = 'AIM CSF 2D';
num_lines = 1;

sName=char(Patient_Name); % assign experimenter values to experimnet
meanLUM=str2num(char('0.5')); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
SFRangecDeg=str2num(char('1 32')); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
paramEsts=str2num(char('3 1.5 0.25 10 0.2')); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
cellSizeDeg=str2num(char('5')); % diameterof cell rings % changed from 2 to 6

nRows=str2num(char('4')); % $ rows and columns on each chart
nCols=str2num(char('4'));
nTrials=number_of_trials; % # charts to run

[width, height]=Screen('DisplaySize', whichScreen); % adaptive screen width measurement in mm
scrnWidthCm = width/10; % converting mm to cm

viewDistCm=str2num(char(view_distance)); % observer's distance from screen
correction=correction; %with (1) or without (2) spectacle correction
testEye=eye; %1=left,2=right, 3=both
feedback=str2num(char('0')); %1=left,2=right, 3=both

interCellGapDeg=1; %
ringLineDeg=0.2;
nextSizeDeg=4; 
ranErrorDeg=45; % random error variable for AIM acuity


%% save data
if correction ==1
    correction='wCorr';
elseif correction ==0
    correction='woCorr';
end

if testEye==1
    testEye='OS'
    txteye='Please cover your right eye and tap to continue';
    soundFile = 'female_cover_right.mp3';
elseif testEye==2
    testEye='OD'
    txteye='Please cover your left eye and tap to continue';
    soundFile = 'female_cover_left.mp3';
elseif testEye==3
    testEye='OU'
    txteye='Look at the screen with both eyes open.Tap to continue';
    soundFile = 'female_use_both_eyes.mp3';
end

testSName=[sName,'_' testEye '_' correction '_' datestr(now,'mm-dd-yyyy HH-MM-SS')];
dataFile=fullfile(cd,'Patient Data', sprintf('%sAIM_CSF2D.mat', testSName));  % file path to unique Mat files
%% screen and stimulus parameters

rng('default'); % seed random number generator
gammaVal=2.22; % gamma on this system
LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background

% feedbackArcRGB=[LMean*1.2 LMean*1.2 LMean*1.2]; % Grayish ring (120 % times LMean) 
% feedbackArcRGB=round(255*(feedbackArcRGB./255).^(1/gammaVal)); % Gamma-corrected target
% 
% gapArcRGB=[LMean*0.5 LMean*0.5 LMean*0.5]; % 50 % of LMean
% gapArcRGB=round(255*(gapArcRGB./255).^(1/gammaVal)); % Gamma-corrected target



feedbackArcRGB=[LMean*0.5 LMean*0.5 LMean*0.5]; % Grayish ring (120 % times LMean)
feedbackArcRGB=round(255*(feedbackArcRGB./255).^(1/gammaVal)); % Gamma-corrected target

gapArcRGB=[LMean*0.1 LMean*0.1 LMean*0.1]; % 50 % of LMean
gapArcRGB=round(255*(gapArcRGB./255).^(1/gammaVal)); % Gamma-corrected target


nAFC=8; % angle for feedback error scoring
saveImage=1; % save stimuli or not

dev = min(GetTouchDeviceIndices([], 1));
if ~isempty(dev)
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
end

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    Screen('Preference','VisualDebugLevel', 0);  % Remove pysch toolbox warning (Red screen)
    %     [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean); %, [0 0 2800 1800]);
    
    if ismac
        rval = kPsychNeedRetinaResolution;
        [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean(1),[],[],[],[],[],rval);%,[0 0 3839, 2160]   ,[0 0 2800 1800]);
    else
        [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean(1));%,[0 0 3839, 2160]   ,[0 0 2800 1800]);
    end
    
    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % height of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    fbRadius=0.75*cellDiamPix/2; % radius of feedback bars
    responseRingPenWidth=pixPerDeg*ringLineDeg; % line width for rresponse ring - TODO justify 0.25 deg
    responseRingRect=[0 0 cellDiamPix-responseRingPenWidth/2 cellDiamPix-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg); % gap betwen cells to minimise interference
    targSizePix=round(cellDiamPix*0.5); % size of targets within cells % changed from 0.8 to 0.5
    srcRect=[0 0 targSizePix targSizePix]; % rect for gabors
    respGapWidthDeg=180/(5*pi); % gap is always 22.9 deg (1/5 of the line width)
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % 'next' bottom, halfway down the screen, away from the left
    
    if ~isempty(dev)
        TouchQueueCreate(windowPtr, dev);
        TouchQueueStart(dev);
        KbReleaseWait;
        HideCursor(windowPtr);
    else
        ShowCursor('Arrow'); % show arrow response cursor
    end
    

    Screen('TextFont',windowPtr,'Arial'); % put some instructions on screen
    Screen('TextSize',windowPtr,100);

    DrawFormattedText(windowPtr, txteye, 'center', 'center', 0, LMean(1)); % put this message on screen
    Screen('Flip', windowPtr); % flip to the information screen
    instructionCommand = audioread(soundFile); % sound file when mouse is clicked
    sound(instructionCommand, 22000)
    WaitSecs(5)

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
    instruction = 'Each cells contain fainter versions of stripes \n Tap on the direction of the stripes in each cell';
    [~, ~, textBox] = DrawFormattedText(windowPtr,  instruction, 'center', 'center',0, LMean(1)); % give instructions to subject
    Screen('Flip', windowPtr); % flip to the information screen

    command = audioread('female_AIMCSFdirection.mp3'); % sound file when mouse is clicked
    sound(command, 22000);
    WaitSecs(7);
    
    Screen('TextSize',windowPtr,48);
    
    nCells=nRows*nCols;
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
    trialSeed=zeros(nTrials);
    
    for trialNo=1:nTrials % set up data structure for each chart
        trialRecord(trialNo) = struct('trialSeed',0,'targSF',[],'targLocs',[],'targContrast', [],'targCSensLevel', [],...
            'targOri', [],'matchOri', [],'oriErr', [],'respCorrect', [],'stimSeen', [],...
            'chartTime', 0,'logsfPeakEsts', 0,'logcsPeakEsts', 0,'bWidthEsts', 0,'minErrEsts', 0 ,'slopeEsts', 0 ); % start empty array
    end
    
    sfRangeLog=linspace(log10(SFRangecDeg(1)), log10(SFRangecDeg(2))); % range of 100 sf values between 1 cpi and nyquist
    csRangeLog=linspace(log10(1/0.001), log10(1/1)); % range of 100 contrast sensitivity values between 0.1% and 100% (inverted b
    [logSF, logCS]=meshgrid(sfRangeLog, csRangeLog); % 2D range of contrasts and possible SFs on this system
    
    xStimCenter=screenCenter(1)+xFactor*(cellDiamPix+interCellGapPix); % x center of each cell on threen for this size target
    yStimCenter=screenCenter(2)+yFactor*(cellDiamPix+interCellGapPix); % y center of each cell on threen for this size target
    
    guessOriRate=45; % for grating, 90 for rings
    
    %% run experiment
    t=1;%counter for speech command
    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials
        
        if trialNo>1 % there has been at least 1 trial
            testSF=[]; % set blank list of all SFs tested
            testCS=[]; % set blank list of all CS levels tested
            oriErr=[]; % blank list of all orientation error levels for all stimuli tested
            for trialSoFar=1:nTrials % work through all trials
                testSF=[testSF trialRecord(trialSoFar).targSF]; % concatenate all SF levels across trials
                testCS=[testCS trialRecord(trialSoFar).targCSensLevel]; % concatenate all CS levels  across trials
                oriErr=[oriErr trialRecord(trialSoFar).oriErr]; % concatenate all Ori Error responses across trials
            end
            fitObj=fitCSFOriErrorSurface_AIMCSF(log10(testSF), testCS, abs(oriErr)); % fit 2D CSF to data

            trialRecord(trialNo).logsfPeakEsts=fitObj.xPeak; % store updated parameter estimates
            trialRecord(trialNo).logcsPeakEsts=fitObj.yPeak;
            trialRecord(trialNo).bWidthEsts=fitObj.bWidth;
            trialRecord(trialNo).minErrEsts=fitObj.minErr;
            trialRecord(trialNo).slopeEsts=fitObj.slope;
        else
            trialRecord(trialNo).logsfPeakEsts=log10(paramEsts(1)); % start with experimenter's parameter estimates
            trialRecord(trialNo).logcsPeakEsts=paramEsts(2);
            trialRecord(trialNo).bWidthEsts=paramEsts(3);
            trialRecord(trialNo).minErrEsts=paramEsts(4);
            trialRecord(trialNo).slopeEsts=paramEsts(5);
        end
        % current estimate of CSF surface with current parameters for the defined range of SFs and CS
        currentCSSurf=CSFOriErrorSurface( logSF, logCS, trialRecord(trialNo).logsfPeakEsts, trialRecord(trialNo).logcsPeakEsts, trialRecord(trialNo).bWidthEsts, trialRecord(trialNo).minErrEsts, trialRecord(trialNo).slopeEsts, guessOriRate, 0);
        
        threshOri= trialRecord(trialNo).minErrEsts+(guessOriRate-trialRecord(trialNo).minErrEsts)/2; % current criterion threshold
        [~, CSLowSFIn]=min(abs(currentCSSurf(:,1)-threshOri)); % index of low SF threshold
        [~, peakSFIn]=min(abs(sfRangeLog-trialRecord(trialNo).logsfPeakEsts)); % index of peak SF
        minOri=atan2d(CSLowSFIn-100,-peakSFIn); % Orientation of left-most threshold
        testCSFAngles=linspace(minOri, 0, nCells); % Orientations in CSF space to test
        
        minThreshOri=max(threshOri-2*trialRecord(trialNo).slopeEsts*threshOri, 0); % orientation error levels to test at each SF/CS combination -  minimum 0
        maxThreshOri=min(threshOri+2*trialRecord(trialNo).slopeEsts*threshOri, guessOriRate); %  maximum random level
        
        testOriStds=linspace(minThreshOri, maxThreshOri, nCells);
        
        [sf, cs]=meshgrid((1:100)-peakSFIn, 99:-1:0); % orientations in CSFs for this chart with current peak SF
        csfOris=atan2d(-cs, sf);
        
        Screen('Flip', windowPtr); % clear screen
        SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
        trialRecord(trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
        rng(trialRecord(trialNo).trialSeed); % seed the random number generator
        
        ranSFOrder=Shuffle(1:nCells); % randomize the SF in each cell
        ranOriOrder=Shuffle(1:nCells); % separately randomize the contrast/ori error level in each cell
        
        trialRecord(trialNo).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect
        
        for cellNum=1:nCells % draw all random orientation rings
            closestOris=abs(csfOris - testCSFAngles(ranSFOrder(cellNum)));
            closestError=abs(currentCSSurf - testOriStds(ranSFOrder(cellNum)));
            closestCombined=closestOris.*closestError;
            [~,idx] = min(closestCombined(:)); % find the minimum value on each row of matrix A
            [rNum,cNum] = ind2sub(size(csfOris),idx);
            trialRecord(trialNo).targSF(cellNum)=10.^sfRangeLog(cNum); % logSF for current target
            trialRecord(trialNo).targCSensLevel(cellNum)=csRangeLog(rNum); % log sensitivity for current target
            trialRecord(trialNo).targContrast(cellNum)=1./10.^csRangeLog(rNum);  % linear contrast for current targer
            testLambda=pixPerDeg/trialRecord(trialNo).targSF(cellNum); % linear wavelength for current target
            
%             myTarg=Gabor2D_AIMCSF(testLambda, 0,360*rand(),cellDiamPix/6,targSizePix,targSizePix/2);
%             myTarg=Gabor2D_AIMCSF(testLambda, 0,360*rand(),cellDiamPix/10,targSizePix,targSizePix/2); %changed (cellDiamPix/10) to remove sharp edges from texture 

            myTarg=Gabor2D_AIMCSF(testLambda, 0,360*rand(),targSizePix/5 ,targSizePix,targSizePix/2); %used targetSizePixels to remove sharp edges from texture

            myTarg=meanLUM+meanLUM*myTarg*trialRecord(trialNo).targContrast(cellNum); % scale to mean luminance - contrast * C
            myTarg=255*myTarg.^(1/gammaVal); % gamma correct
            myTargFloor=floor(myTarg); % base LUT value for each pixel
            myTargResidual=myTarg-myTargFloor; % residual floating point
            MyTargBit=binornd(1, myTargResidual)+myTargFloor; % binomial random sample with probability from residual
            MyTargBit(MyTargBit>255)=255; % catch overspills
            MyTargBit(MyTargBit<0)=0; % catch underspills
            myTex=Screen('MakeTexture', windowPtr,MyTargBit);%, [],[],2); % write image to texture, scaled
            
            testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
            trialRecord(trialNo).targOri(cellNum)=testAngle; % update record for test angle
            
            destRect=CenterRectOnPoint(srcRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
            Screen('DrawTexture', windowPtr, myTex, srcRect,destRect,testAngle); % draw texture to on screen location
            Screen('Close', myTex);
            destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
            Screen('FrameOval', windowPtr,[LMean LMean 255],destRect, responseRingPenWidth, responseRingPenWidth); % draw blue line around cells
        end
        
        Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
        chartStart=tic; % time for each chart
        
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
                                
                                if mx > nextRect(1) && sum(trialRecord(trialNo).stimSeen)==nCells% observer clicked finished word, only when sum of indications made equals nCells
                                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                                    [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                                    mouseOriWRTChosenBox=90+atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                                    clickedExit=1;
                                else
                                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                                    [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                                    mouseOriWRTChosenBox=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                                    destRect=CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback
                                    Screen('FrameOval', windowPtr, feedbackArcRGB, destRect, responseRingPenWidth, responseRingPenWidth); % draw greyish line around cells
                                    Screen('FrameArc', windowPtr,gapArcRGB,destRect,90+mouseOriWRTChosenBox-respGapWidthDeg/2,respGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                                    
                                    beepData = audioread('pop.wav'); % sound file when mouse is clicked
                                    sound(beepData, 50000);
                                    
                                    Screen('FrameArc',windowPtr,gapArcRGB,destRect,-90+mouseOriWRTChosenBox-respGapWidthDeg/2,respGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                                    Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                                end
                                
                                trialRecord(trialNo).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
                                respErr(1)=diff(unwrap([trialRecord(trialNo).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % error of 0 and 180 deg because grating flips at 180
                                respErr(2)=diff(unwrap([trialRecord(trialNo).targOri(respNum),mouseOriWRTChosenBox+180]/180*pi)*180/pi);
                                [~,minErrLoc]=min(abs(respErr)); % observer's error is closer of 2
                                trialRecord(trialNo).oriErr(respNum)=respErr(minErrLoc); % calculate difference betwen actual and reported orientation
                                trialRecord(trialNo).stimSeen(respNum)=1; % mark this cell as response made
                                %computing 'next' buttom
                                if sum( trialRecord(trialNo).stimSeen)==nCells && nTrials~=trialNo %change color of next button to green once all cells have been clicked, say NEXT
                                    Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                                    Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                                    DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [], [],[],[],[],[],nextRect);
                                    
                                    beepData = audioread('good.wav'); % sound file when mouse is clicked
                                    sound(beepData, 48000);
                                    
                                    Screen('TextBackgroundColor', windowPtr, LMean(1)); % change text background colour to default screen color
                                elseif  sum( trialRecord(trialNo).stimSeen)==nCells && nTrials==trialNo %change color of next button to green once all cells have been clicked, say END
                                    Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                                    Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'END' to green
                                    DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [], [],[],[],[],[],nextRect);
                                    
                                    beepData = audioread('good.wav'); % sound file when mouse is clicked
                                    sound(beepData, 48000);
                                    
                                    Screen('TextBackgroundColor', windowPtr, LMean(1)); % change text background colour to default screen color
                                end
                                if abs(trialRecord(trialNo).oriErr(respNum))<(180/(nAFC*2)) % if within nAFC error -
                                    trialRecord(trialNo).respCorrect(respNum)=1; % score correct
                                else
                                    trialRecord(trialNo).respCorrect(respNum)=0; % score incorrect
                                end
                                Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                            else
                                % Below threshold: Kill the blob:
                                blobcol{id} = [];
                            end
                            trialRecord(trialNo).chartTime=toc(chartStart);
                            
                            if feedback==1 %1== feedback is turned ON
                                for cellNum=1:nCells % give feedback
                                    if trialRecord(trialNo).respCorrect(cellNum)==1 % correct response
                                        fbCol= [0 255 0]; % green line
                                    else % incorrect response
                                        fbCol= [255 0 0]; % red line
                                    end
                                    Screen('DrawLine', windowPtr, fbCol, xStimCenter(cellNum)-cosd(trialRecord(trialNo).targOri(cellNum))*fbRadius, yStimCenter(cellNum)-sind(trialRecord(trialNo).targOri(cellNum))*fbRadius, xStimCenter(cellNum)+cosd(trialRecord(trialNo).targOri(cellNum))*fbRadius, yStimCenter(cellNum)+sind(trialRecord(trialNo).targOri(cellNum))*fbRadius ,4);
                                end
                                [~, ~, textBox] = DrawFormattedText(windowPtr, 'Click when ready for next trial:','center', 100, 0); % give instructions to subject
                            end
                            Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
                            
                            if saveImage==0 % saving demo image?
                                myImfileName=sprintf('AIM_CSF%d.jpg', trialNo); % create filename
                                myImage=Screen('GetImage', windowPtr); % grab screnshot
                                imwrite(myImage,myImfileName); % write screenshot to image file
                            end
                        end
                    end
                end
            end
        end

         %% Were the responses significantly different from random chances? If so continue, otherwise show the following screen
            temp_y_90(1:length(trialRecord(trialNo).oriErr(:)))=ranErrorDeg;
            [h,p]=  ttest(abs(trialRecord(trialNo).oriErr(:)),temp_y_90')

            if h  == 0 %  minimum response error is not significantly different from random error either 45 or 90 deg
                % Note: random responses, not different from 90/45
                if(trialNo< nTrials)
                    instcommand1 = audioread('female_AIMCSF_warning.mp3'); % sound file when mouse is clicked
                    sound(instcommand1, 22000)
                    WaitSecs(4)
                end
                trialRecord(trialNo).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response

            elseif h  == 1 %minimum response error is significantly different from random error either 45 or 90 deg
                % min error is significantly different but we are not sure
                % if its larger than 45/90 deg therefore we do t test
                [h1,p1]=  ttest(abs(trialRecord(trialNo).oriErr(:)),ranErrorDeg, 'Tail','right')
                if h1 == 0 %  min error data is not signifantly larger than 90/45 degree based onin right sided tailed t test
                    %Note: regular, valid trial
                    if trialNo==1 && trialNo~=nTrials % if the trial is greater than the 1st trial and smaller than the last trial, it will show the Good job screen
                        % Outputs
                        instcommand3 = audioread('female_welldone_new.mp3'); % sound file when mouse is clicked
                        sound(instcommand3, 22000)
                        WaitSecs(4)                       
                    elseif trialNo >= 2 && trialNo < (nTrials-1)
                        instcommand3 = audioread('female_youaredoinggrest_new.mp3'); % sound file when mouse is clicked
                        sound(instcommand3, 22000)
                        WaitSecs(4)
                        %                     Speak(feedbackTxt);
                    elseif trialNo == (nTrials-1)
                        instcommand4 = audioread('female_onemorechartleft_new.mp3'); % sound file when mouse is clicked
                        sound(instcommand4, 22000)
                        WaitSecs(4)
                    end
                    trialRecord(trialNo).randomchartresp=0; %count number of random responses for each chart %0 == valid response
                elseif h1 == 1  %  min error data is  signifantly larger than 90/45 degree based onin right sided tailed t test
                    %Note: invalid trial
                    if(trialNo< nTrials)
                        instcommand1 = audioread('female_AIMCSF_warning.mp3'); % sound file when mouse is clicked
                        sound(instcommand1, 22000)
                        WaitSecs(4)
                    end
                    trialRecord(trialNo).randomchartresp=1; %count number of random responses for each chart %1 == invalid response beacuse of random response
                end
            end
        
        
    end % end trial loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    if ~isempty(dev)
        TouchQueueStop(dev);
        TouchQueueRelease(dev);
        ShowCursor(windowPtr);
    end
    
    Screen('CloseAll'); % close all windows
    

    save(dataFile,'trialRecord','expDuration','meanLUM','SFRangecDeg','cellSizeDeg','paramEsts', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect', 'examtime'); % save parameters
    output_graph(1) = figure('visible','off');
    testSF=[]; % set blank list of all levels for this SF
    testCS=[]; % set blank list of all levels for this SF
    oriErr=[]; % blank list of all stimuli seen for this SF
    for trialSoFar=1:nTrials % work through all trials
        testSF=[testSF trialRecord(trialSoFar).targSF]; % concatenate all levels
        testCS=[testCS trialRecord(trialSoFar).targCSensLevel]; % concatenate all levels
        oriErr=[oriErr trialRecord(trialSoFar).oriErr]; % concatenate all Ori Error responses
    end
    fitObj=fitCSFOriErrorSurface_AIMCSF(log10(testSF)', testCS', abs(oriErr)'); % fit 2D CSF to data
    currentCSSurf=CSFOriErrorSurface( logSF, logCS, fitObj.xPeak, fitObj.yPeak, fitObj.bWidth, fitObj.minErr, fitObj.slope, guessOriRate, 0); % current estimate of CSF surface
    mesh(sfRangeLog, csRangeLog,currentCSSurf);
    hold on
    scatter3(log10(testSF), testCS,abs(oriErr));
    xlabel('Spatial Frequency'); % label the x axis
    ylabel('Contrast Sensitivity'); % label the x axis
    zlabel('Angular Error [\circ]'); % label y axis
    
    csfCurve=(fitObj.yPeak + log10(0.5).*((1./(((10.^fitObj.bWidth).*log10(2))./2)).*(sfRangeLog - fitObj.xPeak)).^2);
    
    csfAcuity=10.^sfRangeLog(find(csfCurve>0, 1, 'last'));
    AULCSF = trapz(sfRangeLog(csfCurve >= 0), csfCurve(csfCurve >= 0));
    ci=confint(fitObj);
   
    title(sprintf('%s| ID: %s| SF_p_e_a_k %.1fc/deg| CS_p_e_a_k %.1f', testEye, Patient_ID, 10.^fitObj.xPeak, fitObj.yPeak));
    subtitle(['Noise= ' num2str(fitObj.minErr,2), '\circ | Slope= ',num2str(fitObj.slope,2)]); %
    
    hold off
    save(dataFile,'fitObj','csfAcuity', 'AULCSF','trialRecord','expDuration','meanLUM','SFRangecDeg','cellSizeDeg','paramEsts', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect', 'examtime'); % save parameters
    
catch exception
    
    if ~isempty(dev)
        TouchQueueRelease(dev);
    end
    disp(exception.message);
    Screen('CloseAll');
    testSName=[sName,'_' testEye '_' correction '_' ];
    dataFile=fullfile(cd,'Patient Data', sprintf('%sAIM_CSF_2D_FailedRun_%s.mat', testSName, datestr(now,'mm-dd-yyyy HH-MM-SS')));  % file path to unique Mat files
    
    save(dataFile,'trialRecord','meanLUM','SFRangecDeg','cellSizeDeg','paramEsts', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect', 'examtime'); % save parameters
end
end

% 

function [fitobject] =fitCSFOriErrorSurface_AIMCSF(SF, CS, OriErrors)

myFitModel = fittype( 'CSFOriErrorSurface( x, y, xPeak, yPeak, bWidth, minErr, slope, 45, 0)',...
    'independent', {'x','y'}, ...
    'coefficients',{'xPeak','yPeak','bWidth','minErr','slope'},...
    'dependent','z');

fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 1.5 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 2.5 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 3.0 0.5 20 0.4]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
end

function [oriErrSurf, logCSF] =CSFOriErrorSurface( logSF, logCS, xPeak, yPeak, bWidth, minErr, slope, guessRate, plotFigures)

% loSFcpd=0.5; % lower value of CS, c/deg
% hiSFcpd=32; % upper value of CS, c/deg
% xPeak=log10(3); % peak SF
% yPeak=log10(1/0.01); % peak sensitivity
% bWidth=log10(2.0); % SF bandwidth in octaves
% minErr=10; % error for best performance
% slope=0.5; % slope of orientation response
% guessRate=90; % 90 for circles, 45 for grating
% plotFigures=1;
% 
% sfRangeLog=linspace(log10(loSFcpd), log10(hiSFcpd)); % range of 100 sf values between 1 cpi and nyquist
% csRangeLog=linspace(log10(1/1), log10(1/0.001)); % range of 100 contrast sensitivity values between 0.1% and 100%
% [logSF, logCS]=meshgrid(sfRangeLog, csRangeLog); % 2D range of contrasts and possible SFs on this system

logCSF=(yPeak + log10(0.5).*((1./(((10.^bWidth).*log10(2))./2)).*(logSF - xPeak)).^2); % log parabola CSF

oriErrSurf=minErr + (guessRate-minErr).*(0.5-0.5.*erf((logCSF-logCS)./(sqrt(2).*slope)));

if plotFigures==1
    mesh(logSF, logCS, oriErrSurf);
    xlabel('log Spatial Frequency');
    ylabel('log Contrast');
    zlabel('AIM Error (deg)');
end
end


function gabor = Gabor2D_AIMCSF(lambda, theta, phase, stdev, imSize, elCentre)
if nargin==5
    elCentre=round(imSize/2);
end
if size(imSize)==1
    nRows=imSize;
    nCols=imSize;
else
    nRows=imSize(2);
    nCols=imSize(1);
end
if size(elCentre)==1
    yCentre=elCentre;
    xCentre=elCentre;
else
    yCentre=elCentre(2);
    xCentre=elCentre(1);
end
if size(stdev)==1
    stdevY=stdev;
    stdevX=stdev;
else
    stdevY=stdev(1);
    stdevX=stdev(2);
end
angle=theta*pi/180; % orientation deg to radians.
phase=phase*pi/180; % phase deg to radians.
sinX=sin(angle)*((2*pi)/lambda);
cosY=cos(angle)*((2*pi)/lambda);
[x,y]=meshgrid(1-yCentre:nRows-yCentre,1-xCentre:nCols-xCentre);
gabor=exp(-(x.^2)/(2*stdevX.^2)-(y.^2)/(2*stdevY.^2)).*sin(sinX*x+cosY*y+phase);


end

