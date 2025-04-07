% AIM Stereo Bar
% Program to measure stereo acuity threshold based on orientation identification of a depth-defined random dot bar stimulus
% Presents sequence of grids containing random orientation depth-defined bars of a range of stereo disparities
% Disparity of stimuli based on threshold estimate -99%CI to + 99%CI
% Observer clicks perceived orientation of each bar and can change response or leave cells empty - scored as incorrect with random orientation error (90 deg)
% Responses classified as 8AFC and fit with cumulative gaussian, mean and standard error defines size range for next chart
% Data also fit with error function
% 2022 PJB

clear all; close all; commandwindow; % close windows and type into command window
timestamp_experiment_start = datetime('now');
% saveImage=0;
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject Initials', 'Background level (0-1)', 'Peak SF (c/deg)', 'Bar Size (HW deg)','Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)', 'Stereomode- 8-R/B, -1-Polariser', "Disparity sign: 1=crossed (pop-up), 2=uncrossed (pop-down)"};
dlg_title = 'AIM Stereo Bar';
%% lum in candela/ ^2m: 20 for background and target; JSk, Nov 2024
num_lines = 1;
def = {'XX', '0.2', '0', '5 1.25','6', '4','4','4', '70', '67', '8', '1'}; % some defaul values
% def = {'XX', '0.2', '0 2 4', '4 1','6.25', '4','4','2', '70', '80'}; % some defaul values
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % assign experimenter values to experimnet
meanLUM=str2num(char(answer(2,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
fPeak=str2num(char(answer(3,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
barSizeDeg=str2num(char(answer(4,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
cellSizeDeg=str2num(char(answer(5,1))); % diameterof cell rings
nRows=str2num(char(answer(6,1))); % $ rows and columns on each chart
nCols=str2num(char(answer(7,1)));
nTrials=str2num(char(answer(8,1))); % # charts to run
scrnWidthCm=str2num(char(answer(9,1))); % screen width (cm)
viewDistCm=str2num(char(answer(10,1))); % observer's distance from screen
stereoMode= str2num(char(answer(11,1)));
disparitysign= str2num(char(answer(12,1)));

rng('default'); % seed random number generator
gammaVal=2.3; % gamma on this system
LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background
nAFC=8; % for nAFC scoring or angle error
ranErrorDeg=45;
saveImage=0; %1; % save stimuli or not
interCellGapDeg=2; % gap between cells in degree

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    %     [windowPtr, winRect] = PsychImaging('OpenWindow', whichScreen, LMean);

    if stereoMode<0
        crossTalkFactor=-0.1; %-0.095; % StereoCrosstalkReduction - AK edit
        crossLumCorrect=-0.03;
        PsychImaging('AddTask', 'General', 'InterleavedLineStereo', 0);
        LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background
        [windowPtr, winRect] = PsychImaging('OpenWindow', whichScreen, LMean);%, [],  [], [], stereoMode);    % open experimental window with required parameters
        aspectRatio=0.5; % correct for doubling image size
    else
        LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background
        [windowPtr, winRect] = PsychImaging('OpenWindow', whichScreen,  LMean, [0 0 3840 2159],  [], [], stereoMode);    % open experimental window with required parameters
        aspectRatio=1.0; % correct for doubling image size
        crossLumCorrect=0;
    end
    PsychColorCorrection('SetEncodingGamma', windowPtr, 1./gammaVal);
    [sourceFactorOld, destinationFactorOld]=Screen('BlendFunction',windowPtr, GL_SRC_ALPHA, GL_ONE);

    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;

    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    cellRadiusPix=cellDiamPix/2;
    responseRingPenWidth=pixPerDeg*0.25; % line width for rresponse ring - TODO justify 0.25 deg
     %JSk update Nov, 2024
    responseRingRect=[0 0 cellDiamPix+(responseRingPenWidth*3) cellDiamPix+(responseRingPenWidth*3)*aspectRatio]; 
    % responseRingRect=[0 0 cellDiamPix+responseRingPenWidth (cellDiamPix+responseRingPenWidth)*aspectRatio]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);

    barSizeDegH= barSizeDeg(1);
    barSizeDegV= barSizeDeg(2)*aspectRatio;
    barSizeDeg= [barSizeDegH barSizeDegV];
    barSizePix=barSizeDeg*pixPerDeg; % size of bar (H * W) in pixels on this system
    minDisparityRange=0.5; % minium range o contrasts on screen
    nextSizeDeg=cellSizeDeg/2; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg*aspectRatio],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left
    testGapWidthDeg=360/(10*pi); % gap is always 22.9 deg (1/5 of the line width)

    nDots=200; % TODO  justify paramter
    barCoords=[-barSizePix(1)/2 -barSizePix(1)/2 barSizePix(1)/2  barSizePix(1)/2;...
        -barSizePix(2)/2 barSizePix(2)/2 barSizePix(2)/2  -barSizePix(2)/2]; % corners of rectangle for oriented bar

    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,48);                               % put some instructions on screen
    textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nClick on the orientation of each C \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate);
    DrawFormattedText(windowPtr, textToObserver, 100, 100*aspectRatio, 0, LMean); % put this message on screen
    ShowCursor('Arrow'); % show arrow response cursor
    Screen('Flip', windowPtr); % flip to the information screen
    [mx,my,buttons] = GetMouse; % wait for mouse button release before starting experiment
    while any(buttons) % if already down, wait for release
        [mx,my,buttons] = GetMouse(windowPtr);
    end
    while ~any(buttons) % wait for new press
        [mx,my,buttons] = GetMouse(windowPtr);
    end

    nConds=length(fPeak); % how many letter contrasts to test
    lambdaDeg=pixPerDeg./fPeak; % wavelength required for each SF onm this system
    stdPix=lambdaDeg/4.2991; % relationship between center std and wavelength where surround std = 1.6*center std
    destRects=zeros(nDots, 4);

    slopeEsts=2*ones(size(barSizeDeg)); % start with a rough estimate of slope
    minErrEsts=2*ones(size(barSizeDeg)); % start with a rough estimate of internal error
    startDisparityBnds=[log10(1/360) log10(barSizeDeg(2))]'; % 1 arcmin to 1 barwidth

    nCells=nRows*nCols;
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
    trialSeed=zeros(nTrials, nConds);

    for trialNo=1:nTrials % set up data structure for each chart
        for condNo=1:nConds % and for each contrast
            trialRecord(trialNo,condNo) = struct('trialSeed',0,'targDisparityPerLoc',[],'targLocs',[],'targDisparity', [],'targOri', [],'matchOri', [],'oriErr', [],'respCorrect', [],'stimSeen', [],'chartTime', []); % start empty array
        end
    end

    %% run experiment
    % if disparitysign==1
    %         Speak('Please wear the red-blue glasses. Then, look at each cell that contains dots within a ring. Each cell contains a 3D bar that popps out to varying degrees. Use the ring to indicate via click the orientation of each 3D bar.')
    % elseif disparitysign==0
    %         Speak('Please wear the red-blue glasses. Then, look at each cell that contains dots within a ring. Each cell contains a 3D bar that appears behind to varying degrees. Use the ring to indicate via click the orientation of each 3D bar.')
    % end
    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials
        for condNo=1:nConds % work though the list of test contrasts
            chartStart=tic; % time for each chart
            mouseX=[]; mouseY=[]; mouseclick=[];%empty variables for mouse data, JSk edited 2024
            if trialNo>1 % there has been at least 1 chart, fit all data for next chart
                angErr=[]; % set blank list of all levels for this contrast
                testLevel=[]; % set blank list of all levels for this contrast
                stimSeen=[]; % blank list of all stimuli seen for this SF
                respCorrect=[]; % blank list of all stimuli seen for this SF
                for trialSoFar=1:nTrials % work through all trials
                    testLevel=[testLevel trialRecord(trialSoFar,condNo).targDisparityPerLoc]; % concatenate all test sizes
                    angErr=[angErr trialRecord(trialSoFar,condNo).oriErr]; % concatenate all orientation errors
                    stimSeen=[stimSeen trialRecord(trialSoFar,condNo).stimSeen]; % concatenate all responses made
                    respCorrect=[respCorrect  trialRecord(trialSoFar,condNo).respCorrect]; % concatenate all orientation errors
                end
                angErr(stimSeen==0)=45; % assign error on ignored cells to 45 deg
                %                 fitObj=fitPFuncOriErr(testLevel',abs(angErr)',1); % fit with current orientation errors
                fitObj=fitPFuncnAFC_stereo(testLevel, respCorrect, nAFC, 0); % use matlab's fit function to fit cumulative Gaussian to the data
                fitBnds = confint(fitObj, 0.95); % 95% confidence intervals on contrast threshold estimate
                if any(isnan(fitBnds)) % fit failed to estimate acuity and 95% CIs
                    disparityBound(1)=fitObj.thresh-2*fitObj.slope; % lower bound is threshold - 2 std
                    disparityBound(2)=fitObj.thresh+2*fitObj.slope; % upper bound is threshold + 2 std
                else
                    disparityBound(1)=fitBnds(1);
                    disparityBound(2)=fitBnds(2);
                end
                if any(isnan(disparityBound)) % fit failed to estimate slope
                    disparityBound=startDisparityBnds; % use defaults
                end
            else
                disparityBound=startDisparityBnds; % first trial: use defaults
            end
            % set up random dot positions in in each cell
            xDotCenters=cellDiamPix*rand([nCells,nDots])-cellRadiusPix; % random (x,y) positions for all dots in all cells
            yDotCenters=(cellDiamPix*rand([nCells,nDots])-cellRadiusPix); %*aspectRatio;
            dotOri=atan2(yDotCenters, xDotCenters); % orientation of dot wrt center
            xShiftWrap=cos(dotOri+pi)*cellDiamPix; % x and y shifts to wrap dots to other side of center
            yShiftWrap=(sin(dotOri+pi)*cellDiamPix); %*aspectRatio';
            dotEccents=sqrt(xDotCenters.^2+yDotCenters.^2);
            xDotCenters(dotEccents>cellRadiusPix)=xDotCenters(dotEccents>cellRadiusPix)+xShiftWrap(dotEccents>cellRadiusPix); % wrap diametrically dots outside circle
            yDotCenters(dotEccents>cellRadiusPix)=yDotCenters(dotEccents>cellRadiusPix)+yShiftWrap(dotEccents>cellRadiusPix);
            yDotCenters=yDotCenters*aspectRatio;

            minTestLevel=max([startDisparityBnds(1) disparityBound(1)]); % lower 95% CI disparity, minimum 1 arcmin
            maxTestLevel=min([startDisparityBnds(2) disparityBound(2)]); % uppper 95% CI in arcmin and maximum 1 deg
            disparityRange=maxTestLevel-minTestLevel;
            if disparityRange<minDisparityRange % less than 0.5 log range, increase accordingly
                minTestLevel=minTestLevel-(minDisparityRange-disparityRange)/2;
                maxTestLevel=maxTestLevel+(minDisparityRange-disparityRange)/2;
            end

            if fPeak(condNo)>0
                dotSizePix=round(stdPix(condNo)*1.6*8);
                %                 dotSizePix=round(stdPix(condNo)*1.6*8); % size of dog for this condition
                srcRect=[0 0 dotSizePix dotSizePix];
                myDog=fspecial('gaussian',dotSizePix, stdPix(condNo))-fspecial('gaussian',dotSizePix, stdPix(condNo)*1.6);
                myDog=myDog/(max(myDog(:))); % force max 1
                targTex=Screen('MakeTexture', windowPtr,myDog, [],[],2); % convert log image into texture
            else
                dotSizePix=9;
                %                 dotSizePix=10;
                myTarg=makeDisk(dotSizePix,dotSizePix/2);
                srcRect=[0 0 dotSizePix dotSizePix];
                targTex=Screen('MakeTexture', windowPtr,myTarg, [],[],2); % convert log image into texture
            end

            Screen('Flip', windowPtr); % clear screen
            Screen('BlendFunction',windowPtr, GL_SRC_ALPHA, GL_ONE);

            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(trialNo, condNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(trialNo, condNo).trialSeed); % seed the random number generator

            trialRecord(trialNo,condNo).targLocs=Shuffle(1:nCells); % pick random target locations for target contrasts
            trialRecord(trialNo,condNo).targDisparity=linspace(minTestLevel, maxTestLevel, nCells); % log spaced range of contrasts between upper and lower test contrasts
            trialRecord(trialNo,condNo).targDisparityPerLoc(trialRecord(trialNo,condNo).targLocs)=trialRecord(trialNo,condNo).targDisparity; % fill target locations with corresponding contrast

            trialRecord(trialNo,condNo).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,condNo).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,condNo).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,condNo).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

            xStimCenter=screenCenter(1)+xFactor*(cellDiamPix+interCellGapPix); % x center of each cell on threen for this size target
            yStimCenter=screenCenter(2)+yFactor*(cellDiamPix*aspectRatio+interCellGapPix*aspectRatio); % y center of each cell on threen for this size target

            for cellNum=1:nCells % draw all random orientation rings
                pixelDisparity=10.^trialRecord(trialNo,condNo).targDisparityPerLoc(cellNum)*pixPerDeg; % convert log10(disparity) to pixels
                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord(trialNo,condNo).targOri(cellNum)=testAngle; % update record for test angle

                R = [cosd(-testAngle) sind(-testAngle); -sind(-testAngle) cosd(-testAngle)]; % rotation adjustments
                rotCoords = R*(barCoords); % rotate rectangle coordinates
                rotCoords(:,5)=rotCoords(:,1); % complete polygon

                yCentersThisCell=yDotCenters(cellNum,:)+yStimCenter(cellNum); % x position of left eye dots
                xCentersThisCellOS=xDotCenters(cellNum,:); % x position of left eye dots
                xCentersThisCellOD=xDotCenters(cellNum,:); % x position of right eye dots

                targDots=inpolygon(xDotCenters(cellNum,:),yDotCenters(cellNum,:),rotCoords(1,:),rotCoords(2,:)); % find dots inside bar area

                if disparitysign==1
                    xCentersThisCellOS(targDots)=xCentersThisCellOS(targDots)+pixelDisparity/2; % add disparity to dots inside bar area - half to left eye
                    xCentersThisCellOD(targDots)=xCentersThisCellOD(targDots)-pixelDisparity/2; % add disparity to dots inside bar area - half to right eye
                    dis='crossed';
                elseif disparitysign==2
                    xCentersThisCellOS(targDots)=xCentersThisCellOS(targDots)-pixelDisparity/2; % add disparity to dots inside bar area - half to left eye
                    xCentersThisCellOD(targDots)=xCentersThisCellOD(targDots)+pixelDisparity/2; % add disparity to dots inside bar area - half to right eye
                    dis='uncrossed';
                end

                xCentersThisCellOS=xCentersThisCellOS+xStimCenter(cellNum); % x position of left eye dots
                xCentersThisCellOD=xCentersThisCellOD+xStimCenter(cellNum); % x position of right eye dots

                Screen('SelectStereoDrawBuffer', windowPtr, 0); % draw left eye stimulus
                %                  Screen('BlendFunction',windowPtr, GL_SRC_ALPHA, GL_ONE, [1 0 0 0]); % draw red dots
                for dotNum=1:nDots
                    destRects(dotNum,:)= CenterRectOnPoint(srcRect, xCentersThisCellOS(dotNum),yCentersThisCell(dotNum));
                end
                if stereoMode<0
                    destRects(:,2)=destRects(:,2)+dotSizePix/4;
                    destRects(:,4)=destRects(:,2)+dotSizePix*aspectRatio;
                end
                Screen('DrawTextures', windowPtr, targTex, srcRect, destRects', [],[], 0.5); % draw left DoG

                Screen('SelectStereoDrawBuffer', windowPtr, 1); % draw right eye stimulus
                %                 Screen('BlendFunction',windowPtr, GL_SRC_ALPHA, GL_ONE, [0 0 1 0]); % draw blue
                for dotNum=1:nDots
                    destRects(dotNum,:)= CenterRectOnPoint(srcRect, xCentersThisCellOD(dotNum),yCentersThisCell(dotNum));
                end
                if stereoMode<0
                    destRects(:,2)=destRects(:,2)+dotSizePix/4;
                    destRects(:,4)=destRects(:,2)+dotSizePix*aspectRatio;
                end
                Screen('DrawTextures', windowPtr, targTex, srcRect, destRects', [],[], 0.5); % draw left DoG
            end


            Screen('SelectStereoDrawBuffer', windowPtr, 0);
            %              Screen('BlendFunction',windowPtr, GL_ONE, GL_ZERO,[1 0 0 0]); % draw all colors
            for cellNum=1:nCells % draw all random orientation rings
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('FrameOval', windowPtr,35 ,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
            end

            Screen('SelectStereoDrawBuffer', windowPtr, 1);
            %            Screen('BlendFunction',windowPtr, GL_ONE, GL_ZERO,[0 0 1 0]); % draw all colors
            for cellNum=1:nCells % draw all random orientation rings
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                Screen('FrameOval', windowPtr,35 ,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
            end

            [~, ~, textBox] = DrawFormattedText(windowPtr, 'Click on the orientation of each 3D Bar','left', 100*aspectRatio, 0); % give instructions to subject
            Screen('FillOval',windowPtr,[0 100 00], nextRect);  % draw exit button
            %         Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
            DrawFormattedText(windowPtr, sprintf('%s', 'next'), 'center', 'center', 0, [],[],[],[],[],nextRect); %  char(26)

            Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer


            Screen('BlendFunction',windowPtr, sourceFactorOld, destinationFactorOld);

            ShowCursor('Arrow'); % show arrow response cursor
            clickedExit=0; % reset clicked finished
            while clickedExit==0 % keep checking until subject click exit location
                [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
                mouseX=[mouseX mx];
                mouseY=[mouseY my];
                mouseclick=[mouseclick buttons];%empty variables for mouse data, JSk edited 2024

                while any(buttons) % if already down, wait for release
                    [mx,my,buttons] = GetMouse(windowPtr);
                    mouseX=[mouseX mx];
                    mouseY=[mouseY my];
                    mouseclick=[mouseclick buttons];%empty variables for mouse data, JSk edited 2024
                end
                while ~any(buttons) % wait for new press
                    [mx,my,buttons] = GetMouse(windowPtr);
                    mouseX=[mouseX mx];
                    mouseY=[mouseY my];
                    mouseclick=[mouseclick buttons];%empty variables for mouse data, JSk edited 2024
                end

                %  if mx > nextRect(1) % observer clicked finished word
                if mx > nextRect(1) && sum(trialRecord(trialNo).stimSeen)==nCells % observer clicked finished word and all cells were responded to

                    clickedExit=1;
                else
                    my=my*aspectRatio;
                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                    [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                    mouseOriWRTChosenBox=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                    destRect=CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback
                    Screen('SelectStereoDrawBuffer', windowPtr, 0);
                    Screen('FrameOval', windowPtr,LMean*1.1,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
                    Screen('FrameArc',windowPtr,0,destRect, 90+mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                    Screen('FrameArc',windowPtr,0,destRect,-90+mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                    Screen('SelectStereoDrawBuffer', windowPtr, 1);
                    Screen('FrameOval', windowPtr,LMean*1.1,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
                    Screen('FrameArc',windowPtr,0,destRect, 90+mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                    Screen('FrameArc',windowPtr,0,destRect,-90+mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                    Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                end
                trialRecord(trialNo,condNo).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
                respErr(1)=diff(unwrap([trialRecord(trialNo,condNo).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % error of 0 and 180 deg because grating flips at 180
                respErr(2)=diff(unwrap([trialRecord(trialNo,condNo).targOri(respNum),mouseOriWRTChosenBox+180]/180*pi)*180/pi);
                [~,minErrLoc]=min(abs(respErr)); % observer's error is closer of 2
                trialRecord(trialNo,condNo).oriErr(respNum)=respErr(minErrLoc); % calculate difference betwen actual and reported orientation
                %                 trialRecord(trialNo,sizeNum).oriErr(respNum)=diff(unwrap([trialRecord(trialNo,sizeNum).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % calculate difference betwen actual and reported orientation
                trialRecord(trialNo,condNo).stimSeen(respNum)=1; % mark this cell as response made
                %edited by JSk,2024
                %collect mouse data
                trialRecord( trialNo,condNo).MouseX=mouseX; %x coordinates for mouse
                trialRecord( trialNo,condNo).MouseY=mouseY; %y coordinates for mouse
                trialRecord( trialNo,condNo).Mousebuttons=mouseclick; %x coordinates for mouse

                if abs(trialRecord(trialNo,condNo).oriErr(respNum))<(360/(nAFC*2)) % if within nAFC error
                    trialRecord(trialNo,condNo).respCorrect(respNum)=1; % score correct
                else
                    trialRecord(trialNo,condNo).respCorrect(respNum)=0; % score incorrect
                end


            end
            %computing 'next' buttom gets green
            if sum( trialRecord(trialNo,condNo).stimSeen)==nCells && nTrials==trialNo && condNo==nConds%Create END button once all cells have been clicked for all trials and all colors
                Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [], [],[],[],[],[],nextRect); %  char(26)
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green

            elseif sum( trialRecord(trialNo,condNo).stimSeen)==nCells %Create NEXT button once all cells have been clicked for a trial
                Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw NEXT button, green
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [], [],[],[],[],[],nextRect); %  char(26)
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
            end

            trialRecord(trialNo,condNo).chartTime=toc(chartStart);
            if saveImage==1 % saving demo image?
                myImfileName=[ 'AIM_Depth_' num2str(trialNo) '_' num2str(condNo) '_' sName '.jpg']; % create filename
                myImage=Screen('GetImage', windowPtr); % grab screnshot
                imwrite(myImage,myImfileName); % write screenshot to image file
            end


            % ask after first trial whether anything was actually seen in stereoscvpic depth
            % JSk 2024, July
            if trialNo==1 && condNo==1
                Screen('Flip', windowPtr); % clear screen

                while any(buttons) % if already down, wait for release
                    [mx,my,buttons] = GetMouse(windowPtr);
                end
                while ~any(buttons) % wait for new press
                    [mx,my,buttons] = GetMouse(windowPtr);
                end
                % mx > nextRect(1) && sum(trialRecord(trialNo).stimSeen)==nCells % observer clicked finished word and all cells were responded to
                nextSizeDeg=cellSizeDeg; % size of next chart button
                nextRect_seen=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg*aspectRatio],(winRect(3)*0.75)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left


                Screen('FillOval',windowPtr,[32 255 64], nextRect_seen);  % draw exit button, green
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                DrawFormattedText(windowPtr, sprintf('%s', 'Seen Depth'), 'center', 'center', [], [],[],[],[],[],nextRect_seen); %  char(26)
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green

                nextRect_flat=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg*aspectRatio],(winRect(3)*0.35)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left

                Screen('FillOval',windowPtr,[32 255 64], nextRect_flat);  % draw exit button, green
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                DrawFormattedText(windowPtr, sprintf('%s', 'Just flat dots'), 'center', 'center', [], [],[],[],[],[],nextRect_flat); %  char(26)
                Screen('TextBackgroundColor', windowPtr, [32 255 64]); % change text background color of 'NEXT' to green
                Speak('Did you actually see bars in 3D depth or did you see dots on a flat surface?' )
                Screen('Flip', windowPtr); % clear screen

                clickedExit=0; % reset clicked finished
                while clickedExit==0 % keep checking until subject click exit location
                    [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response

                    while any(buttons) % if already down, wait for release
                        [mx,my,buttons] = GetMouse(windowPtr);

                    end
                    while ~any(buttons) % wait for new press
                        [mx,my,buttons] = GetMouse(windowPtr);

                    end

                    %  if mx > nextRect(1) % observer clicked finished word
                    if mx > nextRect_seen(1) % observer clicked on Seen Depth
                        trialRecord(trialNo,condNo).seendepth=1; % 1=yes, seen depth, 0 no
                        clickedExit=1;
                    elseif  mx > nextRect_flat(1) && mx < nextRect_seen(1)% observer clicked on flat, no depth seen
                        trialRecord(trialNo,condNo).seendepth=0; % 1=yes, seen depth, 0 no
                        clickedExit=1;
                    end
                    nextSizeDeg=cellSizeDeg/2; % size of next chart button
                    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg*aspectRatio],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left

                end
            end % end of seen vs not seen depth
            Screen('Close', targTex);

        end % end SFs loop
    end % end trials loop

    expDuration=toc(expStart); % stop timer - how long was experiment?
    Screen('CloseAll'); % close all windows
    testSName=[sName '_' dis ];    % modify sName if subject already exists
    n=1;
    while exist([testSName,'_AIM_Depth.mat'],'file') ~= 0 % check if data file exists with this subject ID
        n=n+1;
        testSName=[sName ,'_' num2str(n) '_' dis]; % if so create new name until unique name found
    end
    dataFile=sprintf('%s_AIM_Depth.mat', testSName);  % file path to unique Mat files
    save(dataFile,'timestamp_experiment_start','trialRecord','expDuration','meanLUM','barSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

    sName=[];
    nConds=1;
    nAFC=8; % for nAFC scoring or angle error
    ranErrorDeg=45;
    figure();
    % figure
    % threshSD=zeros(1, nConds); % set up empty results variable
    % for condNo=1:nConds % plot figure with fits for interleaved conditions - attempt to fit all data
    %     subplot(ceil(nConds^0.5),ceil(nConds/ceil(nConds^0.5)),condNo);
    %
    %     testLevel=[]; % set blank list of all levels for this SF
    %     respCorrect=[]; % blank list of all stimuli seen for this SF
    %     for trialSoFar=1:nTrials % work through all trials
    %         testLevel=[testLevel trialRecord(trialSoFar,condNo).targDisparityPerLoc]; % concatenate all test sizes
    %         respCorrect=[respCorrect  trialRecord(trialSoFar,condNo).respCorrect]; % concatenate all orientation errors
    %     end
    %     [csFitobject,~,threshSD(condNo)]=fitPFuncnAFC_stereo(testLevel, respCorrect, nAFC, 1); % use matlab's fit function to fit cumulative Gaussian to the data
    %     fitBnds = confint(csFitobject, 0.95); % 99% confidence intervals on contrast threshold estimate
    %     xlabel('Disparity level (log10(\circ))'); % label the x axis
    %     ylabel('Proportion Correct'); % label y axis
    % %      title(sprintf('%s, threshold %.2f deg (%.2f %.2f)', sName, 10.^threshSD(condNo), 10.^fitBnds(1), 10.^fitBnds(2)));
    %     title(sprintf(' Threshold %.2f deg (%.2f %.2f)',  10.^threshSD(condNo), 10.^fitBnds(1), 10.^fitBnds(2)));
    %     ylim([0 1]); % fix upper and lower bounds
    % end


    angThreshCS=zeros(1, nConds);
    for condNo=1:nConds % plot figure with fits for interleaved conditions - attempt to fit all data
        subplot(ceil(nConds^0.5),ceil(nConds/ceil(nConds^0.5)),condNo);
        testLevel=[]; % set blank list of all levels for this SF
        respErr=[]; % blank list of all stimuli seen for this SF
        stimSeen=[]; % blank list of all stimuli seen for this SF
        for trialSoFar=1:nTrials % work through all trials
            testLevel=[testLevel trialRecord(trialSoFar,condNo).targDisparityPerLoc]; % concatenate all test sizes
            respErr=[respErr trialRecord(trialSoFar,condNo).oriErr]; % concatenate all orientation errors
            stimSeen=[stimSeen trialRecord(trialSoFar,condNo).stimSeen]; % concatenate all responses made
        end
        respErr(stimSeen==0)=ranErrorDeg; % assign error on ignored cells to 90 deg
        [angleFitobject,~,angThreshCS(condNo)]=fitPFuncOriErr_stereo(testLevel',abs(respErr)',ranErrorDeg,1); % fit cumulative Gaussian to the ori error data
        fitBnds = confint(angleFitobject, 0.95); % 95% confidence intervals on contrast threshold estimate
        xlabel('Disparity level (log10(\circ))'); % label the x axis
        ylabel('Angular Error [\circ]'); % label y axis
        title(sprintf('%s Threshold: %.1f arcsec (%.1f , %.1f)', sName, 10.^angThreshCS(condNo)*3600, 10.^fitBnds(1)*3600, 10.^fitBnds(2)*3600));
        ylim([0 ranErrorDeg*2]); % fix upper and lower bounds
        angErr_mean(condNo)= mean(respErr);%mean error (bias)
        angErr_std(condNo)= std(respErr);%std error (bias)
    end
    save(dataFile,'angErr_mean','angErr_std','timestamp_experiment_start','disparitysign','stereoMode','sName','fPeak','trialRecord','angleFitobject', 'angThreshCS', 'expDuration','meanLUM','barSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

catch exception

    disp(exception.message)
    for i=1:numel(exception.stack)
        exception.stack(i)
    end
    Screen('CloseAll');
    % cd('C:\Users\bexlab\Dropbox\Jan Skerswetat\FInD Things, GA, Faces-repeatibility'); %store data in one folder
    expDuration=toc;

    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIMDepthFailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sAIMDepthFailedRun.mat', testSName);       % file path to Mat files
    % save(dataFile,'trialRecord','testOrder', 'conditionOrder','ranTestOrder','expDuration','faceSize','nRows','nCols','nGenerations','scrnWidthCm','viewDistCm','frameRate','winRect');
    save(dataFile,'exception','timestamp_experiment_start','disparitysign','stereoMode','sName','fPeak','trialRecord', 'expDuration','meanLUM','barSizeDeg','cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

end

% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/10/2022
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject ] = fitPFuncnAFC_stereo( level, response, nAFC, graphData )
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
if graphData
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
    title(sprintf('Thresh. semiconst.=%.4f, 95%% CI (%.4f,%.4f)', thresholdEst, ci(1,1), ci(2,1)));
end
end


% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/10/2022
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject, lettAcuityEquiv] = fitPFuncOriErr_stereo( tSize, mErr,ranErr, graphData,eye,sName, criterion )
% Fit data from nAFC using cumulative gaussian psychometric function.
% Function returns fitted parameters, goodness of fit, and 95%
% confidence intervals for the fitted parameters.
% criterion= which statistical criterion to find the threshold was used.
% Default 1 means half way between noise and max error, 2= 8AFC (Fractstandard), 3=  =ETDRS
% % example data
% level=[0.500000000000000;0.445625469066873;0.397164117362141;0.353972892192069;0.315478672240097;0.281170662595175;0.250593616813636;0.223341796075482;0.199053585276749;0.177406694616788;0.158113883008419;0.140919146563223;0.125594321575479;0.111936056928417;0.0997631157484441;0.125594321575479;0.111936056928417;0.0997631157484441;0.0889139705019463;0.0792446596230558;0.0997631157484441;0.0889139705019463;0.0792446596230558;0.0706268772311378;0.0629462705897085;0.0561009227150983;0.0500000000000001;0.0445625469066874;0.0397164117362142;0.0353972892192070;0.0445625469066874;0.0397164117362142;0.0353972892192070;0.0315478672240097;0.0397164117362142;0.0500000000000001;0.0629462705897085;0.0561009227150983;0.0706268772311379;0.0629462705897085;0.0792446596230559;0.0997631157484443;0.0889139705019464;0.0792446596230559;0.0706268772311380;0.0629462705897086;0.0561009227150984;0.0706268772311380;0.0889139705019465;0.111936056928417];
% response=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;0;1;1;1;1;1;1;1;1;1;0;1;1;1;0;0;0;1;0;1;0;0;1;1;1;1;1;0;0;0;0];
% nAFC=2;
% graphData=1;
% criterion= which threshold criterion was used e.g. 2,4,8,AFC ?

if nargin<4
    eye='';%default etdrs error
end
if nargin<5
    sName='XYZ';%default etdrs error
end
if nargin<6
    criterion=72;%default etdrs error
end
ft = fittype( @(thresh, slope, minErr, x)(minErr + (ranErr-minErr).*(0.5-0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
[fitobject, gof] = fit(tSize, mErr,ft, 'StartPoint', [mean(tSize) std(tSize) 10],  'Upper',[max(tSize) range(tSize) ranErr/2],  'Lower',[min(tSize) 0.05 0]); %Startpoint mean std and arbitary estimate of 10 degree err; max= can't get higher than max size,
% [fitobject, gof] = fit(tSize, mErr,ft);%, 'StartPoint', [mean(tLevel) std(tLevel)]);
thresholdEst=fitobject.thresh;
minErrEst=fitobject.minErr;
ci = confint(fitobject);
expVals=fitobject(tSize); % evaluate fit at test levels

fitSlope=fitobject.slope;

% criterion=72;  % ori error for given # letters per line
%resolve fx for lettAcuityEquiv
% criterion=minErrEst + (90-minErrEst).*(0.5-0.5.*erf((lettAcuityEquiv-thresholdEst)./(sqrt(2).*fitSlope))); % ang error func
% (criterion-minErrEst )/((90-minErrEst))= 0.5-0.5.*erf((lettAcuityEquiv-thresholdEst)./(sqrt(2).*fitSlope)); % ang error func
% ((criterion-minErrEst )/((90-minErrEst))-0.5)/(-0.5) = erf((lettAcuityEquiv-thresholdEst)./(sqrt(2).*fitSlope)); % ang error func
% erfinv(((criterion-minErrEst )/((90-minErrEst))-0.5)/(-0.5)) = (lettAcuityEquiv-thresholdEst)./(sqrt(2).*fitSlope); % ang error func
lettAcuityEquiv=(sqrt(2).*fitSlope)*(erfinv(((criterion-minErrEst )/((ranErr-minErrEst))-0.5)/-0.5))+thresholdEst; % ang error func
chiSqFitReject=chi2gof(sum(((expVals-mErr).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tSize), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)
if graphData==1
    tOriSort=sort(tSize);
    ci95 = predint(fitobject,tOriSort, 0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    scatter(tSize, mErr);
    hold on
    ylim([0,ranErr*2]);
    plot(fitobject, 'r-');
    if ~any(isnan(ci95(:))) % no NaNs in confidence intervals
        plot(tOriSort,ci95,'g--'); % plot 95% confidence interval on fit line
    end
    hold off
    legend HIDE
    box off
    %     set(gca, 'XScale', 'log');
    %title([sprintf('Thresh=%.2f (%.2f,%.2f)', thresholdEst, ci(1,1), ci(2,1))],'FontSize',18);
    %         title([sprintf('  %s| ID:%s | Threshold=  %1.2f ', eye,sName,thresholdEst)],'FontSize',18);
    subtitle(['Noise= ' num2str(fitobject.minErr,2), '\circ  | Slope= ',num2str(fitSlope,2)]); %
    % xlabel('Test Level');
    %ylabel('Angular Error [/circ])','FontSize',18); % label y axis
end
end

% This software is patented and owned by Northeastern University, Boston,
% USA; Patent number: INV-22060
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2022
%% Version 02/10/2022

function myDisk = makeDisk(imSize, circRadius)
% Make image of Standard Landolt C stimulus
% linewidth is 1/5 optotype size
if nargin==1
    circRadius=imSize/2;
end
[X,Y]=meshgrid(linspace(-imSize/2, imSize/2, imSize),linspace(-imSize/2, imSize/2, imSize)); % 2D matrix of radial distances from centre
radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
myDisk=zeros(imSize); % set up background
myDisk(radDist<=circRadius)=1; % make disk of 1's of required size
end

%respIn = [0 0 1 0 1 1 1 0 1 0 1 1 0 1 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 0];
%dataIn =[0.0417 0.0525 0.0661 0.0590 0.0661 0.0624 0.0590 0.0557 0.0624 0.0607 0.0643 0.0624 0.0607 0.0643 0.0624 0.0607 0.0590 0.0573 0.0607 0.0590 0.0573 0.0557 0.0541 0.0573 0.0557 0.0541 0.0525 0.0557 0.0590 0.0624 0.0607];

function arrangedData=arrangeData(dataIn, respIn)
% arangeData
% sort a stream of test levels and responses into psychometric functions
% arrangedData(:,1)= unique testing level
% arrangedData(:,2)= # trials at this level level
% arrangedData(:,3)= # correct/Yes responses at this level
% arrangedData(:,4)= proportion correct/Yes at this level
% arrangedData(:,5)= expected value after rule of succession correction (i.e. add pseudosamples).
% arrangedData(:,6)= beta distribution lower bound (based on Expected value after Rule of Succession correction)
% arrangedData(:,7)= beta distribution upper bound (based on Expected value after Rule of Succession correction)
% testing level

uniqLevels=unique(dataIn);
arrangedData=zeros(length(uniqLevels), 8);
arrangedData(:,1)=uniqLevels';

for index=1:length(uniqLevels)
    trialsThisLevel=(dataIn==uniqLevels(index));            % trials at this level
    arrangedData(index,2)=sum(trialsThisLevel);             % total # trials at this level
    arrangedData(index,3)=sum(trialsThisLevel.*respIn);     % total # correct at this level
end

arrangedData(:,4)=arrangedData(:,3)./arrangedData(:,2);     % proportion correct at each level

% Beta distribution confidence limits:
% expected values:
arrangedData(:,5) = (arrangedData(:,3)+1) ./ (arrangedData(:,2)+2);
nSuccessesE = arrangedData(:,5).*(arrangedData(:,2)+2);
nFailsE = (arrangedData(:,2)+2)-nSuccessesE;
arrangedData(:,6) = betainv(.025,nSuccessesE,nFailsE);
arrangedData(:,7) = betainv(.975,nSuccessesE,nFailsE);

arrangedData(:,8)=((arrangedData(:,4).*(1-arrangedData(:,4))./arrangedData(:,2)).^0.5)+(1./((2*arrangedData(:,2)))); % binomial standard deviation for each lest level, with correction for small sample size
% arrangedData(:,8)=(arrangedData(:,4).*(1-arrangedData(:,4)).*arrangedData(:,2)).^0.5; % binomial standard deviation for each lest level

end





