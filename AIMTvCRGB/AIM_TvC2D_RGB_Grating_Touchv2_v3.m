% AIM TvC Grating
% Program to measure threshold versus contrast function based on
% orientation identification of rotated grating stimulus presented in
% isotropic bandpass filtered noise
% Presents sequence of grids containing random orientation grating stimuli
% of a range of contrasts, on noise pedstals spanning range sepecifed by user
% Contrast increments based on threshold estimate -2 std dev to + 2 stdevs
% Observer clicks perceived orientation of each grating and can change response or leave cells empty - scored as random orientation error (45 deg)
% Response error fit with TvC function and cumulative gaussian psychometric function,
% defines increment contrasts range for next chart
%updated by Jsk January 2026, included Color directions, touchpad, updated
%modl parameters. updated indication rng col, updated size (2 degree
%stimulus size

% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021


clear all; close all; commandwindow; % close windows and type into command window
ExpDate = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm');
% Suppress Psychtoolbox intro and warning screens
Screen('Preference', 'VisualDebugLevel', 0); % Disable splash screen and visual debug warnings
Screen('Preference', 'SkipSyncTests', 1); % Skip sync test warnings
Screen('Preference', 'SuppressAllWarnings', 1); % Suppress additional warnings
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject Initials', 'SF (c/deg)','Pedestal Range (%)', 'Cell Size (deg)',...
    '# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)','Hidden Cell 0 No, 1 yes','With(1) of Without(0) Correction',...
    'Which eye? left=1,right=2,both=3', 'LMS type: 1-BW 2-RG 3-BY 4-L 5-M'};
dlg_title = 'TvC 2D';
num_lines = 1;
def = {'XX', '2', '0.001 0.9', '4', '4','4','5', '28', '40','1','1','3','2'}; % some default values
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % participant ID
targSFcDeg=str2num(char(answer(2,1))); % spatial frequency og grating and noise pedestal
pedRange=str2num(char(answer(3,1))); % contrast range for pedestals
% paramEst=str2num(char(answer(4,1))); % start parameters of dipper function, Boyton et al, 1999 Vision Research
cellSizeDeg=str2num(char(answer(4,1))); % diameter of cell rings; cell size 4 and stimulus size 2 degree, in agreement with anamoloscope stimulus size
nRows=str2num(char(answer(5,1))); % # rows and columns on each chart
nCols=str2num(char(answer(6,1))); % # columns
nTrials=str2num(char(answer(7,1))); % # charts to run
scrnWidthCm=str2num(char(answer(8,1))); % screen width (cm)
viewDistCm=str2num(char(answer(9,1))); % observer's distance from screen
hiddenCell=str2num(char(answer(10,1))); % observer's distance from screen
correction=str2num(char(answer(11,1)));
testEye=str2num(char(answer(12,1))); %1=left,2=right, 3=both
LMSType=str2num(char(answer(13,1))); % LMS type: 1-BW 2-RG 3-BY 4-LIncrLDecr 5-MIncrMDecr

meanLUM=127; % background RGB LUT (will be gamma corrected)
rng('default'); % seed random number generator
gammaVal=2.3; % gamma on this system
LMean=255*(meanLUM/255).^(1/gammaVal); % Gamma-corrected background level
rgbRatios=[1,-0.3425217,0.0061483; 0.1551372, -0.1878286, 1; 1,-0.1423,-0.0163;-1,0.4623,-0.0196]; %[-1,1];
% 2L-M, S, L, M for Surface Pad, vector proportion [0.1677, 0.8599 0.2003 0.2397]

saveImage=0; % save stimuli or not
interCellGapDeg=1; % gap between cells in degrees

if LMSType==1
    LMSType_name='BW';
    % threshEst=[1 .5 .5 .8 2 3 6 8 16];
    % ringCol=[200 200 255];
    % markCol=[255 200 200];
    paramEst=[0.01 0.5];
    ringCol=[LMean*0.9 LMean*0.9 LMean*0.9];
    markCol=[LMean*1.1 LMean*1.1  LMean*1.1];
elseif LMSType==2
    LMSType_name='RG';
    paramEst=[0.05 0.5];
    % threshEst=[3 3 3 5 5 5 9 9 16];
    ringCol=[LMean*0.9 LMean*0.9 LMean*0.9];
    markCol=[LMean*1.1 LMean*1.1  LMean*1.1];
    %     msg="Choose condition";opts=["No glasses" "EnChroma1(Cx1,red)" "EnChroma2(Cx1,blue)" 'NeutralDensityFiltered'];
    %     choice=menu(msg,opts);% 2nd GUI for RG
    % [choice,tf]=listdlg('ListString',{'No glasses','EnChroma1(Cx1,red)','EnChroma2(Cx1,blue)','NeutralDensityFiltered'});
elseif LMSType==3
    LMSType_name='BY';
    paramEst=[0.05 0.5];
    % threshEst=[4 4 4 8 8 8 16 16 16];
    ringCol=[LMean*0.9 LMean*0.9 LMean*0.9];
    markCol=[LMean*1.1 LMean*1.1  LMean*1.1];
elseif LMSType==4
    LMSType_name='L';
    paramEst=[0.01 0.5];
    % threshEst=[4 4 4 8 8 8 16 16 16];
    ringCol=[LMean*0.9 LMean*0.9 LMean*0.9];
    markCol=[LMean*1.1 LMean*1.1  LMean*1.1];
elseif LMSType==5
    LMSType_name='M';
    paramEst=[0.01 0.5];
    % threshEst=[4 4 4 8 8 8 16 16 16];
    ringCol=[LMean*0.9 LMean*0.9 LMean*0.9];
    markCol=[LMean*1.1 LMean*1.1  LMean*1.1];
end

%% save data
if correction ==1
    correction='wCorr';
elseif correction ==0
    correction='woCorr';
end
if testEye==1
    testEye='OS';
    txteye='Please cover your right eye.';
elseif testEye==2
    testEye='OD';
    txteye='Please cover your left eye.';
elseif testEye==3
    testEye='OU';
    txteye='Please use both eyes.';
end

n=1;
testSName=[sName,'_' testEye '_' LMSType_name '_' correction '_' num2str(n)];
% Initialize touchpad device info
dev = min(GetTouchDeviceIndices([], 1));
if ~isempty(dev)
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
end

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    % [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean, [0 0 2800 1800]); % debugging on single screen
    [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean);

    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % height of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    responseRingPenWidth=pixPerDeg*0.25; % line width for rresponse ring
    %JSk update Nov, 2024
    % responseRingRect=[0 0 cellDiamPix+(responseRingPenWidth*3) cellDiamPix+(responseRingPenWidth*3)];

    responseRingRect=[0 0 cellDiamPix-responseRingPenWidth/2 cellDiamPix-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg); % gap between cells in pixels

    targSizePix=round(cellDiamPix*0.5); % target must be smalle than cell to allow rotation without overlap
    testGapWidthDeg=10; % angle for response arc in deg
    testLambda=pixPerDeg./targSFcDeg; % wavelength of stimuli in pixels on this system
    srcRect=[0 0 targSizePix targSizePix]; % rect for stimulus
    mySFFilter=makeFilter(4, targSizePix/testLambda, 0.5, 1, targSizePix); % create bandpass filter at required spatial frequency
    % mySpatialWindow=Gaussian2D(cellDiamPix/6, targSizePix/2, targSizePix);
    % mySpatialWindow=cosWindow2D(targSizePix*0.5, pixPerDeg/2, targSizePix/2,targSizePix/2, targSizePix, targSizePix); % use raised cosine to avoid contrast gradient cosWindow2D(radius, wLength, centreX, centreY, sizeX, sizeY)
    mySpatialWindow=cosWindow2D(targSizePix*0.45, pixPerDeg/2, targSizePix/2,targSizePix/2, targSizePix, targSizePix); % use raised cosine to avoid contrast gradient cosWindow2D(radius, wLength, centreX, centreY, sizeX, sizeY)
    %% updated by JSk, 02/06/2025, 0.45 cuts away the edge
    nextSize=cellDiamPix/2; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSize nextSize],winRect(3)-nextSize,screenCenter(2)); % halfway down the screen, on the right for next chart button
    rgbGabor=zeros(targSizePix,targSizePix,3);
    % Initialize touch or mouse input
    if ~isempty(dev)
        TouchQueueCreate(windowPtr, dev);
        TouchQueueStart(dev);
        KbReleaseWait;
        HideCursor(windowPtr);
    else
        ShowCursor('Arrow'); % Show arrow response cursor
    end
    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,150);                               % put some instructions on screen
    % textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nMinimum logMAR stimulus is %3.4f \\nMaximum logMAR stimulus is %3.4f \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate, minLogMAR, maxLogMAR);
    % DrawFormattedText(windowPtr, textToObserver, 100, 100, 0, backLUT); % put this message on screen
    Screen('DrawText', windowPtr, (txteye), screenCenter(1)/2, screenCenter(2), 0,  LMean(1));
    Speak(txteye)
    vbl=Screen('Flip', windowPtr); % show stimulus and note time
    Screen('TextSize',windowPtr,48); % put some instructions on screen, specify font parameters
    Screen('TextFont',windowPtr,'Arial');

    % textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate);
    % DrawFormattedText(windowPtr, textToObserver, 100, 100, 0, LMean); % put this message on screen
    % ShowCursor('Arrow'); % show arrow response cursor
    % Screen('Flip', windowPtr); % flip to the information screen
    % Wait for touch/click to continue
    if ~isempty(dev)
        evt = TouchEventGet(dev, windowPtr);
        while isempty(evt)
            evt = TouchEventGet(dev, windowPtr);
            if ~isempty(evt)
                if ismember(evt.Type, [1, 2, 3, 4])
                    continue;
                end
            end
        end
    else
        [mx, my, buttons] = GetMouse(windowPtr);
        while any(buttons) % If already down, wait for release
            [mx, my, buttons] = GetMouse(windowPtr);
        end
        while ~any(buttons) % Wait for new press to begin experiment
            [mx, my, buttons] = GetMouse(windowPtr);
        end
        while any(buttons) % Wait for release
            [mx, my, buttons] = GetMouse(windowPtr);
        end
    end

    nCells=nRows*nCols; % total # cells
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size

    for trialNo=1:nTrials % set up data structure for each chart
        trialRecord(trialNo) = struct('trialSeed',0,'INEst',0, 'PsiEst', 0, 'slopeEst', 0,'minErrEst', 0, 'sigLevel',[],'targRanPedLocs',[],'pedLevel', [],'conIncLevel', [],'matchOri', [],'oriErr', [],'oriErrSorted', [],'stimSeen', [],'chartTime', -1); % start empty array
    end
    pedLevels=logspace(log10(pedRange(1)), log10(pedRange(2)), nCells); % range of pedestal levels on each chart, log speced between experimenter's inputs
    sigLevels=linspace(-2, 2, nCells); % range of test incremants, relative to threshold - -2 to +2 stds on current slope estimate

    %% run experiment
    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials

        if trialNo==1 % 1st trial
            trialRecord(trialNo).INEst=paramEst(1); % use experimenter's internal noise guess
            trialRecord(trialNo).PsiEst=paramEst(2); % use experimenter's Psi guess
            trialRecord(trialNo).minErrEst=5; % default estimate of minimum orientation error
            trialRecord(trialNo).slopeEst=0.25; % default estimate of slope - start with a high value to sample high and low contrasts
        else % later trials - use last trial result (correct or incorrect) to adjust signal
            pedLevel=[]; % set blank list of all levels for this contrast
            conLevel=[]; % set blank list of all levels for this contrast
            angErr=[]; % set blank list of all levels for this contrast
            for trialSoFar=1:nTrials % work through all trials
                pedLevel=[pedLevel trialRecord(trialSoFar).pedLevel]; % concatenate all test levels
                conLevel=[conLevel trialRecord(trialSoFar).conIncLevel]; % concatenate all test levels
                angErr=[angErr trialRecord(trialSoFar).oriErrSorted]; % concatenate all ori errs levels
            end
            fitObj=fitTvCOriErrorSurface(pedLevel, conLevel, abs(angErr)); % fit 2D TvC to data
            trialRecord(trialNo).INEst=fitObj.IN; % internal noise fit estimate
            trialRecord(trialNo).PsiEst=fitObj.Psi; % Psi fit estimate
            trialRecord(trialNo).minErrEst=fitObj.minErr; % minimum orientation error  fit estimate
            trialRecord(trialNo).slopeEst=fitObj.slope; % slope fit estimate
        end

        trialRecord(trialNo).pedLevel=pedLevels; % same pedestal contrasts each trial TODO may make this adaptive in future
        trialRecord(trialNo).sigLevel=Shuffle(sigLevels); % range of test points on slope, randomized across pedestal levels

        currentTvCFunc=sqrt((trialRecord(trialNo).INEst.^2+pedLevels.^2).*(1+trialRecord(trialNo).PsiEst))-pedLevels; % TvC function with current paarameters
        trialRecord(trialNo).conIncLevel=10.^(log10(currentTvCFunc)+trialRecord(trialNo).sigLevel*trialRecord(trialNo).slopeEst); % contrast increments at each pedestal - span slope, randomly assigned to pedestals

        trialRecord(trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
        rng(trialRecord(trialNo).trialSeed); % seed the random number generator

        trialRecord(trialNo).targRanPedLocs=Shuffle(1:nCells); % randomize pedestal locations

        xStimCenter=screenCenter(1)+xFactor*(cellDiamPix+interCellGapPix); % x center of each cell on threen for this size target
        yStimCenter=screenCenter(2)+yFactor*(cellDiamPix+interCellGapPix); % y center of each cell on threen for this size target

        mouseX=[]; mouseY=[]; mouseclick=[];%empty variables for mouse data, JSk edited 2024

        for cellNum=1:nCells % work through each location
            if LMSType==1
                ranCell=trialRecord(trialNo).targRanPedLocs(cellNum); % find target for each location
                myTarg=trialRecord(trialNo).conIncLevel(ranCell)*Grating2D(testLambda, 0, 360*rand(), targSizePix); % Gabor2D(testLambda, 0,360*rand(),cellDiamPix/6,targSizePix,targSizePix/2); % grating of required contrast, SF, random phase
                %Grating2D(lambda, theta, phase, imSize)
                myNoise=real(ifft2(fft2(randn(targSizePix)).*mySFFilter)); % new band pass filtered noise sample each cell
                myNoise=myNoise-mean(myNoise(:)); % zero mean
                myNoise=trialRecord(trialNo).pedLevel(ranCell)*(myNoise/max(myNoise(:))); % required contrast
                myTarg=0.5+0.5*mySpatialWindow.*(myTarg+myNoise); % scale to mean luminance - contrast * C
                myTarg=255*myTarg.^(1/gammaVal); % gamma correct
                myTargFloor=floor(myTarg); % base LUT value for each pixel
                myTargResidual=myTarg-myTargFloor; % residual floating point
                MyTargBit=binornd(1, myTargResidual)+myTargFloor; % binomial random sample with probability from residual
                MyTargBit(MyTargBit>255)=255; % catch overspills
                MyTargBit(MyTargBit<0)=0; % catch underspills
                myTex=Screen('MakeTexture', windowPtr,MyTargBit);%, [],[],2); % write image to texture, scaled

                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord( trialNo).targOri(cellNum)=testAngle; % update record for test angle

                destRect=CenterRectOnPoint(srcRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('DrawTexture', windowPtr, myTex, srcRect,destRect,testAngle); % draw texture to on screen location
                Screen('Close', myTex); % clear this texture from memory
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('FrameOval', windowPtr,ringCol*1.03,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells

            else
                ranCell=trialRecord(trialNo).targRanPedLocs(cellNum); % find target for each location
                myTarg=trialRecord(trialNo).conIncLevel(ranCell)*Grating2D(testLambda, 0, 360*rand(), targSizePix); % Gabor2D(testLambda, 0,360*rand(),cellDiamPix/6,targSizePix,targSizePix/2); % grating of required contrast, SF, random phase
                %Grating2D(lambda, theta, phase, imSize)
                myNoise=real(ifft2(fft2(randn(targSizePix)).*mySFFilter)); % new band pass filtered noise sample each cell
                myNoise=myNoise-mean(myNoise(:)); % zero mean
                myNoise=trialRecord(trialNo).pedLevel(ranCell)*(myNoise/max(myNoise(:))); % required contrast

                % myTarg=Grating2D(testLambda, 0,360*rand(),targSizePix);
                % myNoise=real(ifft2(fft2(randn(targSizePix)).*mySFFilter)); % band pass filtered noise
                % myNoise=myNoise-mean(myNoise(:)); % zero mean
                % myNoise=myNoise/max(myNoise(:)); % unit contrast
                % targContrast=10.^trialRecord(trialNo,sfNum).targContrastPerLoc(cellNum);

                if LMSType==2 % red/green modulation
                    meanLUM =0.5;
                    rgbGabor(:,:,1)=mySpatialWindow.*(rgbRatios(1,1)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(1,1)); %L-M
                    rgbGabor(:,:,2)=mySpatialWindow.*(rgbRatios(1,2)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(1,2)); %L-M
                    rgbGabor(:,:,3)=mySpatialWindow.*(rgbRatios(1,3)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(1,3)); %L-M
                elseif LMSType==3 % blue yellow modulation
                    rgbGabor(:,:,1)=mySpatialWindow.*(rgbRatios(2,1)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(2,1)); %S
                    rgbGabor(:,:,2)=mySpatialWindow.*(rgbRatios(2,2)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(2,2)); %S
                    rgbGabor(:,:,3)=mySpatialWindow.*(rgbRatios(2,3)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(2,3)); %S
                    meanLUM =0.5;
                elseif LMSType==4 % L+ L- modulation
                    rgbGabor(:,:,1)=mySpatialWindow.*(rgbRatios(3,1)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(3,1)); %L
                    rgbGabor(:,:,2)=mySpatialWindow.*(rgbRatios(3,2)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(3,2)); %L
                    rgbGabor(:,:,3)=mySpatialWindow.*(rgbRatios(3,3)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(3,3)); %L
                    meanLUM =0.5;
                elseif LMSType==5 % M+ M- modulation
                    rgbGabor(:,:,1)=mySpatialWindow.*(rgbRatios(4,1)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(4,1)); %M
                    rgbGabor(:,:,2)=mySpatialWindow.*(rgbRatios(4,2)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(4,2)); %M
                    rgbGabor(:,:,3)=mySpatialWindow.*(rgbRatios(4,3)*myTarg+(trialRecord(trialNo).pedLevel(ranCell))*myNoise*rgbRatios(4,3)); %M
                    meanLUM =0.5;
                end
                myTarg=255*(meanLUM+meanLUM*rgbGabor).^(1/gammaVal); % gamma correct

                % myTarg=0.5+0.5*mySpatialWindow.*(myTarg+myNoise); % scale to mean luminance - contrast * C
                % myTarg=255*myTarg.^(1/gammaVal); % gamma correct
                myTargFloor=floor(myTarg); % base LUT value for each pixel
                myTargResidual=myTarg-myTargFloor; % residual floating point
                MyTargBit=binornd(1, myTargResidual)+myTargFloor; % binomial random sample with probability from residual
                MyTargBit(MyTargBit>255)=255; % catch overspills
                MyTargBit(MyTargBit<0)=0; % catch underspills
                myTex=Screen('MakeTexture', windowPtr,MyTargBit);%, [],[],2); % write image to texture, scaled

                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord( trialNo).targOri(cellNum)=testAngle; % update record for test angle

                destRect=CenterRectOnPoint(srcRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('DrawTexture', windowPtr, myTex, srcRect,destRect,testAngle); % draw texture to on screen location
                Screen('Close', myTex); % clear this texture from memory
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('FrameOval', windowPtr,ringCol*1.03,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
            end
        end

        [~, ~, textBox] = DrawFormattedText(windowPtr, 'Click on the orientation of each grating:','left', 70, 0); % give instructions to subject

        chartStart=tic; % time for each chart
        vbl = Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer

        ShowCursor('Arrow'); % show arrow response cursor
        clickedExit=0; % reset clicked finished
        buttonPressTime = vbl - 1; % CRITICAL: Initialize to allow immediate first click

        while clickedExit==0 % keep checking until subject clicks exit location

            % Handle touch or mouse input
            if ~isempty(dev)
                while TouchEventAvail(dev)
                    evt = TouchEventGet(dev, windowPtr);
                    if ismember(evt.Type, [1, 2]) % Touch down or move
                        mx = evt.MappedX;
                        my = evt.MappedY;
                        mouseX = [mouseX mx];
                        mouseY = [mouseY my];
                        mouseclick = [mouseclick 1];

                        mouseDistFromEachBox = sqrt((mx - xStimCenter).^2 + (my - yStimCenter).^2);
                        [~, respNum] = min(mouseDistFromEachBox(:));

                        if vbl > buttonPressTime + 0.25 % Wait 250 msec between clicks
                            buttonPressTime = vbl;

                            if mx > nextRect(1) && sum(trialRecord(trialNo).stimSeen) == nCells
                                clickedExit = 1;
                                break;
                            else
                                mouseOriWRTChosenBox = atan2d(my - yStimCenter(respNum), mx - xStimCenter(respNum));
                                destRect = CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum));
                                Screen('FrameOval', windowPtr, ringCol, destRect, responseRingPenWidth, responseRingPenWidth);
                                Screen('FrameArc', windowPtr, markCol, destRect, 90 + mouseOriWRTChosenBox - testGapWidthDeg/2, testGapWidthDeg, responseRingPenWidth, responseRingPenWidth);
                                Screen('FrameArc', windowPtr, markCol, destRect, -90 + mouseOriWRTChosenBox - testGapWidthDeg/2, testGapWidthDeg, responseRingPenWidth, responseRingPenWidth);

                                trialRecord(trialNo).matchOri(respNum) = mouseOriWRTChosenBox;
                                respErr(1) = diff(unwrap([trialRecord(trialNo).targOri(respNum), mouseOriWRTChosenBox] / 180 * pi) * 180 / pi);
                                respErr(2) = diff(unwrap([trialRecord(trialNo).targOri(respNum), mouseOriWRTChosenBox + 180] / 180 * pi) * 180 / pi);
                                [~, minErrLoc] = min(abs(respErr));
                                trialRecord(trialNo).oriErr(respNum) = respErr(minErrLoc);
                                trialRecord(trialNo).oriErrSorted(trialRecord(trialNo).targRanPedLocs(respNum)) = trialRecord(trialNo).oriErr(respNum);
                                trialRecord(trialNo).stimSeen(respNum) = 1;

                                % Update NEXT button
                                if sum(trialRecord(trialNo).stimSeen) == nCells && nTrials ~= trialNo
                                    Screen('FillOval', windowPtr, [32 255 64], nextRect);
                                    DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', 255, [],[],[],[],[], nextRect);
                                elseif sum(trialRecord(trialNo).stimSeen) == nCells && nTrials == trialNo
                                    Screen('FillOval', windowPtr, [32 255 64], nextRect);
                                    DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', 255, [],[],[],[],[], nextRect);
                                end

                                vbl = Screen('Flip', windowPtr, [], 1);
                            end
                        end
                    end
                end
            else
                % Mouse input processing
                [mx, my, buttons] = GetMouse(windowPtr);
                mouseX = [mouseX mx];
                mouseY = [mouseY my];
                mouseclick = [mouseclick any(buttons)];

                mouseDistFromEachBox = sqrt((mx - xStimCenter).^2 + (my - yStimCenter).^2);
                [~, respNum] = min(mouseDistFromEachBox(:));

                if any(buttons) && vbl > buttonPressTime + 0.25
                    buttonPressTime = vbl;

                    if mx > nextRect(1) && sum(trialRecord(trialNo).stimSeen) == nCells
                        clickedExit = 1;
                    else
                        mouseOriWRTChosenBox = atan2d(my - yStimCenter(respNum), mx - xStimCenter(respNum));
                        destRect = CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum));
                        Screen('FrameOval', windowPtr, ringCol, destRect, responseRingPenWidth, responseRingPenWidth);
                        Screen('FrameArc', windowPtr, markCol, destRect, 90 + mouseOriWRTChosenBox - testGapWidthDeg/2, testGapWidthDeg, responseRingPenWidth, responseRingPenWidth);
                        Screen('FrameArc', windowPtr, markCol, destRect, -90 + mouseOriWRTChosenBox - testGapWidthDeg/2, testGapWidthDeg, responseRingPenWidth, responseRingPenWidth);

                        trialRecord(trialNo).matchOri(respNum) = mouseOriWRTChosenBox;
                        respErr(1) = diff(unwrap([trialRecord(trialNo).targOri(respNum), mouseOriWRTChosenBox] / 180 * pi) * 180 / pi);
                        respErr(2) = diff(unwrap([trialRecord(trialNo).targOri(respNum), mouseOriWRTChosenBox + 180] / 180 * pi) * 180 / pi);
                        [~, minErrLoc] = min(abs(respErr));
                        trialRecord(trialNo).oriErr(respNum) = respErr(minErrLoc);
                        trialRecord(trialNo).oriErrSorted(trialRecord(trialNo).targRanPedLocs(respNum)) = trialRecord(trialNo).oriErr(respNum);
                        trialRecord(trialNo).stimSeen(respNum) = 1;

                        % Update NEXT button
                        if sum(trialRecord(trialNo).stimSeen) == nCells && nTrials ~= trialNo
                            Screen('FillOval', windowPtr, [32 255 64], nextRect);
                            DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', 255, [],[],[],[],[], nextRect);
                        elseif sum(trialRecord(trialNo).stimSeen) == nCells && nTrials == trialNo
                            Screen('FillOval', windowPtr, [32 255 64], nextRect);
                            DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', 255, [],[],[],[],[], nextRect);
                        end

                        vbl = Screen('Flip', windowPtr, [], 1);
                    end
                end
            end

            % CRITICAL: Update vbl even when no click occurs
            WaitSecs(0.001); % Small delay to prevent CPU overuse
            vbl = GetSecs(); % Update current time

            % Collect mouse/touch data
            trialRecord(trialNo).MouseX = mouseX;
            trialRecord(trialNo).MouseY = mouseY;
            trialRecord(trialNo).Mousebuttons = mouseclick;
        end % End while clickedExit

        trialRecord( trialNo).chartTime=toc(chartStart); % time to complete this chart
        if saveImage==1 % saving demo image?
            myImfileName=sprintf('AIM_TvC_2D%d.jpg', trialNo); % create filename
            myImage=Screen('GetImage', windowPtr); % grab screenshot
            imwrite(myImage,myImfileName); % write screenshot to image file
        end
        vbl=Screen('Flip', windowPtr); % show stimulus and note time

        % Wait for button release
        if ~isempty(dev)
            TouchEventFlush(dev);
            WaitSecs(0.5);
        else
            [mx, my, buttons] = GetMouse(windowPtr);
            while any(buttons)
                [mx, my, buttons] = GetMouse(windowPtr);
            end
        end
    end % End trials loop

    expDuration = toc(expStart);

    if ~isempty(dev)
        TouchQueueStop(dev);
        TouchQueueRelease(dev);
        ShowCursor(windowPtr);
    end

    expDuration=toc(expStart); % stop timer -record total duration of test
    Screen('CloseAll'); % close all windows
    n=1;
    while exist([testSName,'AIM_TvC_2D.mat'],'file') ~= 0 % check if data file exists with this subject ID
        n=n+1;
        testSName=[sName,'_' testEye '_' LMSType_name '_' correction '_' num2str(n)];
    end
    dataFile=sprintf('%sAIM_TvC_2D.mat', testSName);  % file path to unique Mat files
    save(dataFile,'ExpDate','trialRecord','expDuration','meanLUM','targSFcDeg','pedRange','paramEst', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','pixPerDeg','frameRate','winRect'); % save parameters before fit
    pedLevel=[]; % set blank list of all levels for this test
    conLevel=[];
    angErr=[];
    for trialSoFar=1:nTrials % work through all trials
        pedLevel=[pedLevel trialRecord(trialSoFar).pedLevel]; % concatenate all pedestal levels
        conLevel=[conLevel trialRecord(trialSoFar).conIncLevel]; % concatenate all contrast increment levels
        angErr=[angErr trialRecord(trialSoFar).oriErrSorted]; % concatenate all ori errs
    end
    fitObj=fitTvCOriErrorSurface(pedLevel, conLevel, abs(angErr)); % fit 2D TvC to data
    figure(); % open a new figure
    cIncLevels=logspace(log10(0.001), log10((1)), 100); % range of contrast increment levels for visualization surface
    [pedVal, cIncVal]=meshgrid(pedLevels, cIncLevels); % 2D range of test conditions on this system
    TvCOriErrorSurface(pedVal, cIncVal, fitObj.IN, fitObj.Psi, fitObj.slope, fitObj.minErr, 45, 1); % plot surface with final fit estimates
    hold on
    scatter3(pedLevel, conLevel,abs(angErr));% scatter raw data on figure
    set(gca, 'YScale', 'log'); % set axes
    set(gca, 'XScale', 'log');
    xlabel('Pedestal Contrast');
    ylabel('Increment Threshold');
    zlabel('AIM Error [\circ]');
    title(['IN:' num2str(fitObj.IN) ' Psi:' num2str(fitObj.Psi)])
    save(dataFile,'ExpDate','fitObj', 'trialRecord','expDuration','meanLUM','targSFcDeg','pedRange','paramEst', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','pixPerDeg','frameRate','winRect'); % save parameters with fit


catch exception
    Screen('CloseAll');
    testSName=[sName,'_' testEye '_'  LMSType_name '_' correction '_' num2str(n)];
    n=1;
    while exist([testSName,'AIM_TvC' '_2D_FailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,'_' testEye '_'  LMSType_name '_' correction '_' num2str(n)];
    end
    dataFile=sprintf('%sAIM_TvCColor_2D_FailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'exception','ExpDate','trialRecord','expDuration','meanLUM','targSFcDeg','pedRange','paramEst', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','pixPerDeg','frameRate','winRect'); % save parameters before fit
end

function grating = Grating2D(lambda, theta, phase, imSize)

angle=theta*pi/180; % orientation deg to radians.
phase=phase*pi/180; % phase deg to radians.
sinX=sin(angle)*((2*pi)/lambda);
cosY=cos(angle)*((2*pi)/lambda);

if size(imSize)<2
    imSize(2)=imSize(1);
end

[x,y]=meshgrid(1:imSize(1),1:imSize(2));

grating=sin(sinX*x+cosY*y+phase);
end

function madeFilter = makeFilter(filtType, fPeak, bWdth, alpha, filtSize)
% Function to make a range of basic filters
% filtType=5; filtSize=512; bWdth=10*pi/180; fPeak=64; alpha=1;
% update svn.
if length(filtSize)==1
    filtSize(2)=filtSize(1);
end
filtRadius=round(filtSize/2);

[X,Y]=meshgrid(-filtRadius(2):-filtRadius(2)+filtSize(2)-1,-filtRadius(1):-filtRadius(1)+filtSize(1)-1);                      % 2D matrix of radial distances from centre
radDist=(X.^2+Y.^2).^0.5;                                                                                   % radial distance from centre
radDist(filtRadius(1),filtRadius(2))=0.5;                                                                         % avoid log or divide by zero

if      filtType == 1   madeFilter = exp(-((log(2)*(abs(log(radDist/fPeak))).^3)/((bWdth*log(2)).^3)));     % isotropic log exponential
elseif  filtType == 2   madeFilter = radDist.^-alpha;                                                       % isotropic 1/f^-alpha
elseif  filtType == 3   fPeak=fPeak*pi/180; bWdth=bWdth*pi/180;                                             % convert degrees to radians
    angDist=atan2(-Y, X);                                                               % orientation filter - angular dist
    sintheta = sin(angDist);
    costheta = cos(angDist);
    ds = sintheta * cos(fPeak) - costheta * sin(fPeak);                                 % Difference in sine
    dc = costheta * cos(fPeak) + sintheta * sin(fPeak);                                 % Difference in cosine
    dtheta = abs(atan2(ds,dc));                                                         % Absolute angular distance
    madeFilter = exp((-dtheta.^2) / (2 * bWdth^2));                                     % Calculate the angular filter component

    fPeak=fPeak+pi;                                                             % 180 deg offset in +ve TFs
    ds = sintheta * cos(fPeak) - costheta * sin(fPeak);                         % Difference in sine
    dc = costheta * cos(fPeak) + sintheta * sin(fPeak);                         % Difference in cosine
    dtheta = abs(atan2(ds,dc));                                                     % Absolute angular distance
    madeFilter = madeFilter+exp((-dtheta.^2) / (2 * bWdth^2));                                     % Calculate the angular filter component

    % elseif  filtType == 4   radDist=log2(radDist);                                                             % isotropic log Cosine - Nasanen et al 1998
    %                         madeFilter = 0.5*(1+cos(pi*(radDist-log2(fPeak))));
    %                         madeFilter(radDist>(log2(fPeak)+1))=0;
    %                         madeFilter(radDist<=(log2(fPeak)-1))=0;
elseif  filtType == 4   radDist=log2(radDist);                                                             % isotropic log Cosine - Nasanen et al 1998
    madeFilter = 0.5*(1+cos(pi*(radDist-log2(fPeak))));
    madeFilter(radDist>(log2(fPeak)+1))=0;
    madeFilter(radDist<=(log2(fPeak)-1))=0;
elseif  filtType == 5   madeFilter = exp(-((log2(radDist)-log2(fPeak)).^2)/(2*(bWdth))^2);                  % log Gaussian
elseif  filtType == 6   madeFilter = exp(-((radDist-fPeak).^2)/(2*bWdth^2));                                %  Gaussian
else                    if fPeak<0 madeFilter = radDist>=abs(fPeak);                                        % default: hat box, -ve cut off = high-pass,
else madeFilter = radDist<=fPeak;                                                    % +ve cut-off = low pass hat box
end
end
% subplot(2,1,1); imagesc(madeFilter);
% subplot(2,1,2); plot(madeFilter(filtRadius,:))
madeFilter=fftshift(madeFilter);                                                                            % shift FFT for compatability with FFT2 (0 cpi is top left etc.)
madeFilter(1,1)=0;
end

% cosWindow
% cWin=cosWindow(radius, wLength, centreX, centreY, sizeX, sizeY)
%
% Makes an image containing a cosinusoidal edged disk.
% radius - radius of window disk to centre of cosine smoothing function
% wLength - peak to trough distance of cosine smoothing (striclty, half the cosine wavelength)
% SizeX and SizeY - the x and y size of the image

function cWin=cosWindow2D(radius, wLength, centreX, centreY, sizeX, sizeY)

[x,y]=meshgrid(1:sizeX,1:sizeY);
radDist=((x-centreX).^2+(y-centreY).^2).^0.5; % distance from centre
cWin=radDist;
cWin(radDist<=(radius-wLength/2))=1; % middle set to 1
cWin(radDist>=(radius+wLength/2))=0; % outside set to 0
cWin(radDist>(radius-wLength/2) & radDist<(radius+wLength/2))=0.5+0.5*cos(pi/2+pi*(cWin(radDist>(radius-wLength/2) & radDist<(radius+wLength/2))-radius)/wLength);

end


function [fitobject] =fitTvCOriErrorSurface(PedVal, CSVal, OriErrors)

% myFitModel = fittype( 'TvCOriErrorSurface( x, y, IN, Psi, slope, minErr, 45, 0)',...
myFitModel = fittype( '45 - (45-minErr).*(0.5-0.5.*erf((log10(sqrt((IN.^2+x.^2).*(1+Psi))-x)-log10(y))./(sqrt(2).*slope)))',...
    'independent', {'x','y'}, ...
    'coefficients',{'IN','Psi','slope','minErr'},...
    'dependent','z');

fitobject = fit( [PedVal(:), CSVal(:)], OriErrors(:), myFitModel, 'StartPoint', [0.02 0.5 0.1 5], 'Lower',[0.01 0.01 0.01 1], 'Upper',[0.1 2 0.2 20]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);

end

function oriErrSurf =TvCOriErrorSurface( pedVal, cIncVal, IN, Psi, slope, minErr, guessRate, plotFigures)

% pedRange=[.01 .64];
% nCells=16;
% pedVal=logspace(log10(pedRange(1)), log10(pedRange(2)), nCells); % range of pedestal levels
% cIncVal=logspace(log10(0.002), log10(0.32), 100); % range of contrast increment levels
% [pedVal, cIncVal]=meshgrid(pedVal, cIncVal); % 2D range of contrasts and possible SFs on this system
% IN = 0.02;
% Psi = 0.5;
% slope = 0.1;
% minErr = 5;
% guessRate=90; % 90 for circles, 45 for grating
% plotFigures=1;
%
% TvCfunc=sqrt((IN.^2+pedVal.^2).*(1+Psi))-pedVal;

% oriErrSurf=minErr + (guessRate-minErr).*(0.5-0.5.*erf((TvCfunc-cIncVal)./(sqrt(2).*slope)));

% oriErrSurf=minErr + (guessRate-minErr).*(0.5-0.5.*erf((log10(TvCfunc)-log10(cIncVal))./(sqrt(2).*slope)));
oriErrSurf=guessRate - (guessRate-minErr).*(0.5-0.5.*erf((log10(sqrt((IN.^2+pedVal.^2).*(1+Psi))-pedVal)-log10(cIncVal))./(sqrt(2).*slope)));

if plotFigures==1
    mesh(pedVal, cIncVal, oriErrSurf);
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlabel('log(Pedestal Contrast)');
    ylabel('log(Increment Threshold)');
    zlabel('AIM Error (deg)');
end
end
