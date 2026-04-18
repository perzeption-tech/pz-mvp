% AIM CSF
% Program to measure contrast threshold based on orientation identification of rotated Gabor stimulus
% Presents sequence of grids containing random orientation Gabor stimuli of a range of contrasts
% Contrast of stimuli based on threshold estimate -2 std dev to + 2 stdevs
% Observer clicks perceived orientation of each Gabor and can change response or leave cells empty - scored as incorrect with random orientation error (90 deg)
% Responses classified as 8AFC and fit with cumulative gaussian, mean and stnadard error defines size range for next chart
% Data also fit with error function, based on angular difference between true and reported Gabor orientatio
%
% 2022  PJB
% WhiteVal=252.5; % luminance of LUT 255
% BlackVal=1.5; % luminance of LUT 0
% meanLUT=LumDesired/WhiteVal*255;
clear all; close all; commandwindow; % close windows and type into command window

whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject Initials','1-BW 2-RG 3-BY', 'Background level (0-1)', 'SF range (c/deg)','Parameter Estimate (%)','Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)'};
dlg_title = 'FInD Acuity';
num_lines = 1;
def = {'XX', '1', '0.5', '1 32','3 1.5 0.25 10 0.2', '7', '4','4','2', '70', '60'}; % some defaul values
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % assign experimenter values to experimnet
LMSType=str2num(char(answer(2,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
meanLUM=str2num(char(answer(3,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
SFRangecDeg=str2num(char(answer(4,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
paramEsts=str2num(char(answer(5,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
cellSizeDeg=str2num(char(answer(6,1))); % diameterof cell rings
nRows=str2num(char(answer(7,1))); % $ rows and columns on each chart
nCols=str2num(char(answer(8,1)));
nTrials=str2num(char(answer(9,1))); % # charts to run
scrnWidthCm=str2num(char(answer(10,1))); % screen width (cm)
viewDistCm=str2num(char(answer(11,1))); % observer's distance from screen

rng('default'); % seed random number generator
gammaVal=2.0; % gamma on this system
LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background
rgbRatios=[1.0, -0.102998, -0.00317653; -1.0, 0.478886, -0.0109738; 0.122854, -0.197033, 1.0]; % RGB ratios for LG 32"

nAFC=8; % angle for error scoring
saveImage=0; % save stimuli or not
interCellGapDeg=0.75;
SFsupportFactor=1.5;

[y, Fs] = audioread('mixkit-game-ball-tap-2073.wav');


try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    [windowPtr, winRect]=Screen('OpenWindow', whichScreen, LMean);%, [0 0 3000 2000]);
    
    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    fbRadius=0.75*cellDiamPix/2;
    responseRingPenWidth=pixPerDeg*0.25; % line width for rresponse ring - TODO justify 0.25 deg
    responseRingRect=[0 0 cellDiamPix-responseRingPenWidth/2 cellDiamPix-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);
    targSizePix=round(cellDiamPix*0.8);
    respGapWidthDeg=180/(5*pi); % gap is always 22.9 deg (1/5 of the line width)
    srcRect=[0 0 targSizePix targSizePix];
    rgbGabor=zeros(targSizePix,targSizePix,4);
    nextSizeDeg=cellSizeDeg/2; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left
    
    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,48);                               % put some instructions on screen
    textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate);
    DrawFormattedText(windowPtr, textToObserver, 100, 100, 0, LMean); % put this message on screen
    ShowCursor('Arrow'); % show arrow response cursor
    Screen('Flip', windowPtr); % flip to the information screen
    [mx,my,buttons] = GetMouse; % wait for mouse button release before starting experiment
    while any(buttons) % if already down, wait for release
        [mx,my,buttons] = GetMouse(windowPtr);
    end
    while ~any(buttons) % wait for new press
        [mx,my,buttons] = GetMouse(windowPtr);
    end
    
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
            fitObj=fitCSFOriErrorSurface(log10(testSF), testCS, abs(oriErr)); % fit 2D CSF to data
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
        minOri=atan2d(CSLowSFIn-100,-peakSFIn); % Orientation of left most threshold
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
            
            myTarg=Gabor2D(testLambda, 0,360*rand(),cellDiamPix/6,targSizePix,targSizePix/2);
            targContrast=trialRecord(trialNo).targContrast(cellNum);
            if LMSType==1  % luminance defined
                rgbGabor(:,:,1)=targContrast*myTarg; %L-M
                rgbGabor(:,:,2)=targContrast*myTarg; %L-M
                rgbGabor(:,:,3)=targContrast*myTarg; %L-M
            elseif LMSType==2 % red/green modulation
                lmsContrast=sqrt((targContrast.^2)/2);
                rgbGabor(:,:,1)=lmsContrast*(rgbRatios(1,1)-rgbRatios(2,1))*myTarg; %L-M
                rgbGabor(:,:,2)=lmsContrast*(rgbRatios(1,2)-rgbRatios(2,2))*myTarg; %L-M
                rgbGabor(:,:,3)=lmsContrast*(rgbRatios(1,3)-rgbRatios(2,3))*myTarg; %L-M
                % cone contrast = sqrt(L^2 + M^2 + S^2); Kim et al (2018). A Normative Data Set for the Clinical Assessment of Achromatic and Chromatic Contrast Sensitivity Using a qCSF Approach. Investigative Opthalmology & Visual Science, 58(9), 3628.
            else % blue yellow modulation
                % Cc^2 = 3*Clms^2; (where L=M=S=1)
                % (Cc^2)/2 = Clms^2;
                % Clms=sqrt((Cc^2)/2)
                lmsContrast=sqrt((targContrast.^2)/3);
                rgbGabor(:,:,1)=lmsContrast*((rgbRatios(1,1)+rgbRatios(2,1))/2-rgbRatios(3,1))*myTarg; %(L+M)/2-S
                rgbGabor(:,:,2)=lmsContrast*((rgbRatios(1,2)+rgbRatios(2,2))/2-rgbRatios(3,2))*myTarg; %(L+M)/2-S
                rgbGabor(:,:,3)=lmsContrast*((rgbRatios(1,3)+rgbRatios(2,3))/2-rgbRatios(3,3))*myTarg; %(L+M)/2-S
            end
                        
            myTargRGB=meanLUM+meanLUM*rgbGabor; % scale to mean luminance - contrast * C
            myTargRGB=255*myTargRGB.^(1/gammaVal); % gamma correct
            myTargFloor=floor(myTargRGB); % base LUT value for each pixel
            myTargResidual=myTargRGB-myTargFloor; % residual floating point
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
            Screen('FrameOval', windowPtr,200,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
        end
        
        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
        DrawFormattedText(windowPtr, sprintf('%s', 'next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
        
        Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
        chartStart=tic; % time for each chart
        
        ShowCursor('Arrow'); % show arrow response cursor
        clickedExit=0; % reset clicked finished
        while clickedExit==0 % keep checking until subject click exit location
            [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
            while any(buttons) % if already down, wait for release
                [mx,my,buttons] = GetMouse(windowPtr);
            end
            while ~any(buttons) % wait for new press
                [mx,my,buttons] = GetMouse(windowPtr);
            end
            %                 if mx < textBox(3) && my < textBox(4) % observer clicked finished area
            %             if mx < 1.5*cellDiamPix && my < 1.5*cellDiamPix % observer clicked finished area
            if mx > nextRect(1) % observer clicked finished word
                clickedExit=1;
            else
                soundsc(y,Fs);
                mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                mouseOriWRTChosenBox=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback
                Screen('FrameOval', windowPtr,[200 200 255],destRect, responseRingPenWidth, responseRingPenWidth); % draw blue line around cells
                Screen('FrameArc',windowPtr,[255 200 200],destRect,90+mouseOriWRTChosenBox-respGapWidthDeg/2,respGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                Screen('FrameArc',windowPtr,[255 200 200],destRect,-90+mouseOriWRTChosenBox-respGapWidthDeg/2,respGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
            end
            trialRecord(trialNo).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
            respErr(1)=diff(unwrap([trialRecord(trialNo).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % error of 0 and 180 deg
            respErr(2)=diff(unwrap([trialRecord(trialNo).targOri(respNum),mouseOriWRTChosenBox+180]/180*pi)*180/pi);
            [~,minErrLoc]=min(abs(respErr));
            %                 trialRecord(trialNo).oriErr(respNum)=diff(unwrap([trialRecord(trialNo).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % calculate difference betwen actual and reported orientation
            trialRecord(trialNo).oriErr(respNum)=respErr(minErrLoc); % calculate difference betwen actual and reported orientation
            trialRecord(trialNo).stimSeen(respNum)=1; % mark this cell as response made
            if abs(trialRecord(trialNo).oriErr(respNum))<(180/(nAFC*2)) % if within nAFC error -
                trialRecord(trialNo).respCorrect(respNum)=1; % score correct
            else
                trialRecord(trialNo).respCorrect(respNum)=0; % score incorrect
            end
        end
        trialRecord(trialNo).chartTime=toc(chartStart);
        
        %         for cellNum=1:nCells % give feedback
        %             if trialRecord(trialNo).respCorrect(cellNum)==1 % correct response
        %                 fbCol= [0 255 0];
        %             else
        %                 fbCol= [255 0 0];
        %             end
        %             Screen('DrawLine', windowPtr, fbCol, xStimCenter(cellNum)-cosd(trialRecord(trialNo).targOri(cellNum))*fbRadius, yStimCenter(cellNum)-sind(trialRecord(trialNo).targOri(cellNum))*fbRadius, xStimCenter(cellNum)+cosd(trialRecord(trialNo).targOri(cellNum))*fbRadius, yStimCenter(cellNum)+sind(trialRecord(trialNo).targOri(cellNum))*fbRadius ,4);
        %         end
        %         Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
        
        if saveImage % saving demo image?
            myImfileName=sprintf('AIM_CSF%d.jpg', trialNo); % create filename
            myImage=Screen('GetImage', windowPtr); % grab screnshot
            imwrite(myImage,myImfileName); % write screenshot to image file
        end
    end % end SFs loop
    
    % end % end trials loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    Screen('CloseAll'); % close all windows
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIM_CSF_2D.mat'],'file') ~= 0 % check if data file exists with this subject ID
        n=n+1;
        testSName=[sName,num2str(n)]; % if so create new name until unique name found
    end
    dataFile=sprintf('%sAIM_CSF_2D.mat', testSName);  % file path to unique Mat files
    save(dataFile,'trialRecord','expDuration','meanLUM','SFRangecDeg','cellSizeDeg','paramEsts', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters
    
    testSF=[]; % set blank list of all levels for this SF
    testCS=[]; % set blank list of all levels for this SF
    oriErr=[]; % blank list of all stimuli seen for this SF
    for trialSoFar=1:nTrials % work through all trials
        testSF=[testSF trialRecord(trialSoFar).targSF]; % concatenate all levels
        testCS=[testCS trialRecord(trialSoFar).targCSensLevel]; % concatenate all levels
        oriErr=[oriErr trialRecord(trialSoFar).oriErr]; % concatenate all Ori Error responses
    end
    fitObj=fitCSFOriErrorSurface(log10(testSF)', testCS', abs(oriErr)'); % fit 2D CSF to data
    currentCSSurf=CSFOriErrorSurface( logSF, logCS, fitObj.xPeak, fitObj.yPeak, fitObj.bWidth, fitObj.minErr, fitObj.slope, guessOriRate, 0); % current estimate of CSF surface
    mesh(sfRangeLog, csRangeLog,currentCSSurf);
    hold on
    scatter3(log10(testSF), testCS,abs(oriErr));
    xlabel('Spatial Frequency'); % label the x axis
    ylabel('Contrast Sensitivity'); % label the x axis
    zlabel('Angular Error (deg)'); % label y axis
    
    %     [~, csfIdx]=min(abs(currentCSSurf - 0.75*guessOriRate)); % calculate threshold value for each SF
    csfCurve=(fitObj.yPeak + log10(0.5).*((1./(((10.^fitObj.bWidth).*log10(2))./2)).*(sfRangeLog - fitObj.xPeak)).^2);
    %       csfCurve=csRangeLog(csfIdx);
    csfAcuity=10.^sfRangeLog(find(csfCurve>0, 1, 'last'));
    AULCSF = trapz(sfRangeLog(csfCurve >= 0), csfCurve(csfCurve >= 0));
    
    title(sprintf('%s, SFpeak %.1fc/deg CSpeak %.1f Acuity %.1fc/deg AULCSF  %.2f', sName, 10.^fitObj.xPeak, fitObj.yPeak, csfAcuity, AULCSF));
    hold off
    save(dataFile,'fitObj','csfAcuity', 'AULCSF','trialRecord','expDuration','meanLUM','SFRangecDeg','LMSType','cellSizeDeg','paramEsts', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters
    
catch
    Screen('CloseAll');
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIM_CSF_2D_FailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sAIM_CSF_2D_FailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'trialRecord','expDuration','meanLUM','SFRangecDeg','cellSizeDeg','paramEsts','LMSType', 'cellSizeDeg','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters
end



