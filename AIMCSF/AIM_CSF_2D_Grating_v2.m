% AIM CSF
% Program to measure contrast threshold based on orientation identification of rotated Gabor stimulus
% Presents sequence of grids containing random orientation Gabor stimuli of a range of spatial frequencies and contrasts
% Spatial frequecy range M*N equal log spaced between lowest presentable SF % (based on Gaussian stdev) and current estimate of CSF acuity
% Contrast of stimuli based on threshold estimate -2 std dev to + 2 stdevs at each test SF
% Observer clicks perceived orientation of each Gabor and can change response or leave cells empty - scored as incorrect with random orientation error (90 deg)
% Data fit with error function, based on absolute angular difference between true and reported Gabor orientatio
%
% % This software is patented and owned by Northeastern University, Boston,
% % USA; Patent number: PCT/US23/12595
% % Invented and developed by Peter J. Bex and Jan Skerswetat, 2022clear all; close all; commandwindow; % close windows and type into command window
ExpDate = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm')

whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject ID', 'Background level (0-1)', 'Gabor Sigma (deg)','Parameter Estimate (%)','Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)'};
dlg_title = 'AIM CSF 2D';
num_lines = 1;
def = {'XX', '0.5', '1','3 2 0.25 10 0.5', '7', '4','4','3', '70', '60'}; % some defaul values
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % Participant ID
meanLUM=str2num(char(answer(2,1))); % Display mean luminance
gaussSigmaDeg=str2num(char(answer(3,1))); % Range of Spatial frequencies to evaluate
paramEsts=str2num(char(answer(4,1))); % Starting Estimates of Log parabola CSF
cellSizeDeg=str2num(char(answer(5,1))); % diameterof cell rings // Specified in degrees of viewing angle
nRows=str2num(char(answer(6,1))); % $ rows and columns on each chart
nCols=str2num(char(answer(7,1)));
nTrials=str2num(char(answer(8,1))); % # charts to run
scrnWidthCm=str2num(char(answer(9,1))); % screen width (cm)
viewDistCm=str2num(char(answer(10,1))); % observer's distance from screen

rng('default'); % seed random number generator
gammaVal=2.3; % gamma on this system
LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background

nAFC=8; % angle for correct/incorrect error scoring
saveImage=1; % save stimuli or not
interCellGapDeg=0.75;
S 

%[y, Fs] = audioread('G:\My Drive\bexlabfunctions\trunk\AIM\mixkit-game-ball-tap-2073.wav'); % feedback sounds
% [y, Fs] = audioread('C:\Users\pbex\OneDrive\Documents\MATLAB\mixkit-game-ball-tap-2073.wav'); % feedback sounds

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
    responseRingPenWidth=pixPerDeg*0.25; % line width for response ring - TODO justify 0.25 deg
    responseRingRect=[0 0 cellDiamPix-responseRingPenWidth/2 cellDiamPix-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);
    gaussSigmaPix=gaussSigmaDeg*pixPerDeg; % standard deviation of gabor envelope in pixels
    targSizePix=round(gaussSigmaPix*6); % truncate target +/- 3 stdevs
    SFRangecDeg=[pixPerDeg/(2*gaussSigmaPix) pixPerDeg/(2*sqrt(2))]; % min testable SF is 1 cycle per 2 sigma maximum testable SF is wavelength 2*sqrt 2 pixels to allow for 45deg rotated gratings

    respGapWidthDeg=180/(5*pi); % gap is always 22.9 deg (1/5 of the line width matches Landolt C)
    srcRect=[0 0 targSizePix targSizePix];
    nextSizeDeg=cellSizeDeg/2; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left
    % put some instructions on screen
    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,48);
    textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nRange %.1f to %.1f c/deg \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate, SFRangecDeg(1),SFRangecDeg(2));
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
        trialRecord(trialNo) = struct('trialSeed',0,'targlogSF',[],'targLocs',[],'targContrast', [],'targCSensLevel', [],'sigmaLevel', [],...
            'targOri', [],'matchOri', [],'oriErr', [],'oriErrSorted', [],'respCorrect', [],'stimSeen', [],...
            'chartTime', 0,'logsfPeakEsts', 0,'logcsPeakEsts', 0,'bWidthEsts', 0,'minErrEsts', 0 ,'slopeEsts', 0 ); % start empty array
    end

    sfRangeLog=linspace(log10(SFRangecDeg(1)), log10(SFRangecDeg(2))); % range of 100 sf values between 1 cpi and nyquist
    csRangeLog=linspace(log10(1/0.001), log10(1/1)); % range of 100 contrast sensitivity values between 0.1% and 100% (inverted b
    [logSF, logCS]=meshgrid(sfRangeLog, csRangeLog); % 2D range of contrasts and possible SFs on this system

    xStimCenter=screenCenter(1)+xFactor*(cellDiamPix+interCellGapPix); % x center of each cell on threen for this size target
    yStimCenter=screenCenter(2)+yFactor*(cellDiamPix+interCellGapPix); % y center of each cell on threen for this size target

    sigmaLevels=linspace(-2, 2, nCells); % range of test incremants, relative to threshold - -2 to +2 stds on current slope estimate
    guessOriRate=45; % for grating, 90 for rings

    %% run experiment
    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials

        if trialNo>1 % there has been at least 1 trial
            testlogSF=[]; % set blank list of all SFs tested
            testCS=[]; % set blank list of all CS levels tested
            oriErr=[]; % blank list of all orientation error levels for all stimuli tested
            for trialSoFar=1:nTrials % work through all trials
                testlogSF=[testlogSF trialRecord(trialSoFar).targlogSF]; % concatenate all SF levels tested across trials
                testCS=[testCS trialRecord(trialSoFar).targCSensLevel]; % concatenate all CS levels tested across trials
                oriErr=[oriErr trialRecord(trialSoFar).oriErrSorted]; % concatenate all Ori Error responses across trials
            end
            fitObj=fitCSFOriErrorSurface(testlogSF, testCS, abs(oriErr)); % fit 2D CSF to data
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

        Screen('Flip', windowPtr); % clear screen
        SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center

        trialRecord(trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
        rng(trialRecord(trialNo).trialSeed); % seed the random number generator

        currentCSF=(trialRecord(trialNo).logcsPeakEsts + log10(0.5).*((1./(((10.^trialRecord(trialNo).bWidthEsts).*log10(2))./2)).*(sfRangeLog - trialRecord(trialNo).logsfPeakEsts)).^2); % current CSF curve, high res
        upperTestSF=min([SFRangecDeg(2) 10.^sfRangeLog(find(currentCSF>0, 1, 'last'))]); % current CSF acuity estimate, or upper testable SF
        trialRecord(trialNo).targlogSF=linspace(log10(SFRangecDeg(1)), log10(upperTestSF), nCells); % spatial frequencies between lowest SF and acuity
        testLambda=pixPerDeg./(10.^trialRecord(trialNo).targlogSF); % linear wavelength for SF on scurrent chart

        currentCSF=(trialRecord(trialNo).logcsPeakEsts + log10(0.5).*((1./(((10.^trialRecord(trialNo).bWidthEsts).*log10(2))./2)).*(trialRecord(trialNo).targlogSF - trialRecord(trialNo).logsfPeakEsts)).^2); % CSF for current test SFs
        trialRecord(trialNo).sigmaLevel=Shuffle(sigmaLevels); % range of test points on current slope, randomized across SFs
        trialRecord(trialNo).targCSensLevel=currentCSF+trialRecord(trialNo).sigmaLevel*trialRecord(trialNo).slopeEsts; % test contrast sensitivity +/- current slope 
        trialRecord(trialNo).targCSensLevel(trialRecord(trialNo).targCSensLevel<0)=0; % catch contrast greater than 100%
        trialRecord(trialNo).targCSensLevel(trialRecord(trialNo).targCSensLevel>log10(1/0.001))=log10(1/0.001); % catch contrast less than .1%
        trialRecord(trialNo).targContrast=1./(10.^trialRecord(trialNo).targCSensLevel); % linear contrasts per SF on current chart 

        trialRecord(trialNo).targLocs=Shuffle(1:nCells); % randomize SF locations

        trialRecord(trialNo).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

        for cellNum=1:nCells % draw all random orientation rings
            ranCell=trialRecord(trialNo).targLocs(cellNum); % find target for each location
            myTarg=Gabor2D(testLambda(ranCell), 0,360*rand(),gaussSigmaPix,targSizePix,targSizePix/2); % gabor of required wavelength
            myTarg=meanLUM+meanLUM*myTarg*trialRecord(trialNo).targContrast(ranCell); % scale to mean luminance - contrast * C
            myTarg=255*myTarg.^(1/gammaVal); % gamma correct
            myTargFloor=floor(myTarg); % base LUT value for each pixel
            myTargResidual=myTarg-myTargFloor; % residual floating point
            MyTargBit=binornd(1, myTargResidual)+myTargFloor; % binomial random sample with probability from residual
            MyTargBit(MyTargBit>255)=255; % catch overspills
            MyTargBit(MyTargBit<0)=0; % catch underspills
            myTex=Screen('MakeTexture', windowPtr,MyTargBit);%, [],[],2); % write image to texture, scaled

            testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
            trialRecord(trialNo).targOri(cellNum)=testAngle; % update record for test angle

            destRect=CenterRectOnPoint(srcRect, xStimCenter(cellNum), yStimCenter(cellNum)); % location for this grating
            Screen('DrawTexture', windowPtr, myTex, srcRect,destRect,testAngle); % draw texture to on screen location
            Screen('Close', myTex); % release texture
            destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % location for response ring
            Screen('FrameOval', windowPtr,200,destRect, responseRingPenWidth, responseRingPenWidth); % draw response ring around cells
        end

        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
        DrawFormattedText(windowPtr, sprintf('%s', 'next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)

        Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
        chartStart=tic; % start time for each chart

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
            if mx > nextRect(1) % observer clicked finished word
                clickedExit=1;
            else
               % soundsc(y,Fs); % play beep when clicked
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
            [~,minErrLoc]=min(abs(respErr)); % closest response angle to taget angle
            trialRecord(trialNo).oriErr(respNum)=respErr(minErrLoc); % calculate difference betwen actual and reported orientation
            trialRecord(trialNo).oriErrSorted(trialRecord(trialNo).targLocs(respNum))=trialRecord(trialNo).oriErr(respNum); % store response in correct pedestal location
            trialRecord(trialNo).stimSeen(respNum)=1; % mark this cell as response made
            if abs(trialRecord(trialNo).oriErr(respNum))<(180/(nAFC*2)) % if within nAFC error -
                trialRecord(trialNo).respCorrect(respNum)=1; % score correct
            else
                trialRecord(trialNo).respCorrect(respNum)=0; % score incorrect
            end
        end
        trialRecord(trialNo).chartTime=toc(chartStart); % save time to finish this chart

        if saveImage % saving demo image?
            myImfileName=sprintf('AIM_CSF_2D_Chart_%d.jpg', trialNo); % create filename
            myImage=Screen('GetImage', windowPtr); % grab screnshot
            imwrite(myImage,myImfileName); % write screenshot to image file
        end
    end % end SFs loop

    % end % end trials loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    Screen('CloseAll'); % close all windows
    testSName=sName; % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIM_CSF_2D.mat'],'file') ~= 0 % check if data file exists with this subject ID
        n=n+1;
        testSName=[sName,num2str(n)]; % if so create new name until unique name found
    end
    dataFile=sprintf('%sAIM_CSF_2D.mat', testSName);  % file path to unique Mat files
    save(dataFile,'trialRecord','expDuration','meanLUM','gaussSigmaDeg', 'SFRangecDeg','cellSizeDeg','paramEsts','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

    testlogSF=[]; % set blank list of all levels for this SF
    testCS=[]; % set blank list of all levels for this SF
    oriErr=[]; % blank list of all stimuli seen for this SF
    for trialSoFar=1:nTrials % work through all trials
        testlogSF=[testlogSF trialRecord(trialSoFar).targlogSF]; % concatenate all levels
        testCS=[testCS trialRecord(trialSoFar).targCSensLevel]; % concatenate all levels
        oriErr=[oriErr trialRecord(trialSoFar).oriErrSorted]; % concatenate all Ori Error responses
    end
    fitObj=fitCSFOriErrorSurface(testlogSF', testCS', abs(oriErr)'); % fit 2D CSF to data
    currentCSSurf=CSFOriErrorSurface( logSF, logCS, fitObj.xPeak, fitObj.yPeak, fitObj.bWidth, fitObj.minErr, fitObj.slope, guessOriRate, 0); % current estimate of CSF surface
    mesh(sfRangeLog, csRangeLog, currentCSSurf); % plot current CSF surface
    hold on
    scatter3(testlogSF, testCS,abs(oriErr)); % scatter raw data for each cell response
    xlabel('Spatial Frequency'); % label the x axis
    ylabel('Contrast Sensitivity'); % label the x axis
    zlabel('Angular Error (deg)'); % label y axis

    csfCurve=(fitObj.yPeak + log10(0.5).*((1./(((10.^fitObj.bWidth).*log10(2))./2)).*(sfRangeLog - fitObj.xPeak)).^2); % Calculate this CSF curvce at high resolution
    csfAcuity=10.^sfRangeLog(find(csfCurve>0, 1, 'last')); % estimate final CSF acuity
    AULCSF = trapz(sfRangeLog(csfCurve >= 0), csfCurve(csfCurve >= 0)); % Calculate final AULCSF

    title(sprintf('%s, SFpeak %.1fc/deg CSpeak %.1f Acuity %.1fc/deg AULCSF  %.2f', sName, 10.^fitObj.xPeak, fitObj.yPeak, csfAcuity, AULCSF));
    hold off
    save(dataFile,'ExpDate','fitObj','csfAcuity', 'AULCSF','trialRecord','expDuration','meanLUM','gaussSigmaDeg','SFRangecDeg','cellSizeDeg','paramEsts','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

catch exception
    Screen('CloseAll');
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIM_CSF_2D_FailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sAIM_CSF_2D_FailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'ExpDate','exception','trialRecord','expDuration','meanLUM','gaussSigmaDeg', 'SFRangecDeg','cellSizeDeg','paramEsts','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters
end

function [fitobject] =fitCSFOriErrorSurface(SF, CS, OriErrors)

myFitModel = fittype( 'CSFOriErrorSurface( x, y, xPeak, yPeak, bWidth, minErr, slope, 45, 0)',...
     'independent', {'x','y'}, ...
     'coefficients',{'xPeak','yPeak','bWidth','minErr','slope'},...
     'dependent','z');

% fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 1.5 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
% fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 2.6 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.4 0 0.05 0.5 0.1], 'Upper',[1 2.5 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);

 
end





%% CSF 2D function AIM
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
    surf(logSF, logCS, oriErrSurf);
     % hold on 
     % scatter3(logSF, logCS, oriErrSurf);
    xlabel('log Spatial Frequency','FontSize',20)
    ylabel('log Contrast','FontSize',20)
    zlabel('Angular Error [ \circ ]','FontSize',20)
    title('AIM CSF+','FontSize',18);
  %  subtitle('Participant with albinism','FontSize',18);

xlim([0 1.6]);


end
end
% 
% % This software is patented and owned by Northeastern University, Boston,
% % USA; Patent number: INV-22060
% % Invented and developed by Peter J. Bex and Jan Skerswetat, 2022
% %% Version 02/10/2022
% 
% function [fitobject] =fitCSFOriErrorSurface(SF, CS, OriErrors)
% 
% myFitModel = fittype( 'CSFOriErrorSurface( x, y, xPeak, yPeak, bWidth, minErr, slope, 45, 0)',...
%      'independent', {'x','y'}, ...
%      'coefficients',{'xPeak','yPeak','bWidth','minErr','slope'},...
%      'dependent','z');
% 
% % fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 1.5 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
% % fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 2.5 0.5 20 0.5]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
%  fitobject = fit( [SF(:), CS(:)], OriErrors(:), myFitModel, 'StartPoint', [0.5 1.5 0.25 10 0.2], 'Lower',[-0.3 0 0.1 1 0.1], 'Upper',[0 3.0 0.5 20 0.4]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);
% 
% 
% end


function gabor = Gabor2D(lambda, theta, phase, stdev, imSize, elCentre)

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
