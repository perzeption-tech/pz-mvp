
% AIM VA
% Program to measure visual acuity based on orientation identification of rotated Landolt C stimulus
% Presents sequence of grids containing random orientation C stimuli of a  range of sizes
% Size of stimuli based on threshold estimate -2 std dev to + 4 stdevs
% Observer clicks perceived orientation of each C and can change response or leave cells empty - scored as incorrect with random orientation error (90 deg)
% Responses classified as 8AFC and fit with cumulative gaussian, mean and stnadard error defines size range for nxt chart
% Data also fit with error function
% 2022  PJB
clear all; close all; commandwindow; % clear all memory, close windows and type into command window
ExpDate = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm')

whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject Initials', 'Contrast','Kernel width (deg)','Background Luminance (0-1)','Ring Luminance (0-1)', 'Ring Thickness (deg)', 'Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)','With(1) of Without(0) Correction'};%,'Which eye? left=1,right=2,both=3'};
dlg_title = 'AIM VA Blur';
num_lines = 1;
%def = {'XX', '1','0.0 0.1 0.2 0.4 0.8','0.1','0.6', '0.25', '10', '4','4','3', '70', '40','0','3'}; % some default values
def = {'XX', '0.97','0.0','0.5','0.7', '0.1', '1.25', '4','4','4', '70', '400','0'}; % some default values

answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % assign experimenter values to experimnet
tContrast=str2num(char(answer(2,1))); % Weber contrast of letters
blurWidth=str2num(char(answer(3,1)));
meanLUM=str2num(char(answer(4,1))); % Weber contrast of letters

ringLUM=str2num(char(answer(5,1))); % Weber contrast of letters
ringLineDeg=str2num(char(answer(6,1))); % Weber contrast of letters
cSize=str2num(char(answer(7,1))); % diameter of cell rings
nRows=str2num(char(answer(8,1))); % $ rows and columns on each chart
nCols=str2num(char(answer(9,1)));
nTrials=str2num(char(answer(10,1))); % # charts to run
scrnWidthCm=str2num(char(answer(11,1))); % screen width (cm)
viewDistCm=str2num(char(answer(12,1))); % observer's distance from screen
correction=str2num(char(answer(13,1))); %with (1) or without (2) spectacle correction
%testEye=str2num(char(answer(14,1))); %1=left,2=right, 3=both

rng('default'); % seed random number generator
gammaVal=2.3; % gamma on this system
% backLUT=255*meanLUM.^(1/gammaVal); % Gamma-corrected background
%meanLUT=194; %taken from 2024 AIM Acuity HL VA
meanLUT=2; %low lum JSk,09/07/2024

backLUT =255*(meanLUT./255).^(1/gammaVal); % Gamma-corrected background
%backLUT=255*meanLUM; % Gamma-corrected background
targLum=meanLUM-(tContrast.*meanLUM); % luminance of target 0-1
% targLUT=round(255*targLum.^(1/gammaVal)); % Gamma-corrected target
targLUT=round(255*targLum); % Gamma-corrected target
% respLUT=round(255*ringLUM.^(1/gammaVal)); % Gamma-corrected target
respLUT=round(255*ringLUM)*0.1; % Gamma-corrected target
ranNumSeed=1;
nAFC=8; % angle for error scoring
saveImage=1; % save stimuli or not
interCellGapDeg=0.25;

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    %     [windowPtr, winRect]=Screen('OpenWindow', whichScreen, backLUT);

    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible'); % request 32 bit per pixel for high res contrast
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    [windowPtr, winRect]=PsychImaging('OpenWindow', whichScreen, backLUT);
    PsychColorCorrection('SetEncodingGamma', windowPtr, 1/gammaVal);

    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    imSize=round(cSize*pixPerDeg); % size of each cell in pixels on this system
    srcRect=[0 0 imSize imSize]; % rect for this sized cell
    responseRingPenWidth=pixPerDeg*ringLineDeg; % line width for rresponse ring - TODO justify 0.25 deg
    responseRingRect=[0 0 imSize-responseRingPenWidth/2 imSize-responseRingPenWidth/2]; % size of the cell border
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);


    pixSizeArcmin=60/pixPerDeg; % size of eacdh pixel in arcmin
    minLogMAR=log10(pixSizeArcmin); % minimum target size logMAR
    maxLogMAR=log10((pixSizeArcmin*imSize)/12); % maximum target size logMAR
    %    maxLogMAR=log10(0.2*(scrnHeightPix/pixSizeArcmin)); % maximum target size logMAR - letter is full screen height
    nextSize=imSize/2; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSize nextSize],winRect(3)-nextSize/2,screenCenter(2)); % halfway down the screen, away from the left

    % create file to save data - make sure not to overwrite existing file
    storeRNG=rng(round(sum(100*clock)));
    randomeye=randperm(3);%randomise the order of eyes
    rng(ranNumSeed);

    for i=1:length(randomeye)%run through 3 eye conditions
        testEye=randomeye(i);

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

        testSName=[sName,'_' testEye '_' correction '_' ];
        n=1;
        while exist([testSName,'AIM_VA_Blur.mat'],'file') ~= 0 % check if data file exists with this subject ID
            n=n+1;
            testSName=[sName,num2str(n)]; % if so create new name until unique name found
        end
        dataFile=sprintf('%sAIM_VA_Blur.mat', testSName);  % file path to unique Mat files

        Screen('TextFont',windowPtr,'Arial');
        Screen('TextSize',windowPtr,150);                               % put some instructions on screen
        % textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nMinimum logMAR stimulus is %3.4f \\nMaximum logMAR stimulus is %3.4f \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate, minLogMAR, maxLogMAR);
        % DrawFormattedText(windowPtr, textToObserver, 100, 100, 0, backLUT); % put this message on screen
        Screen('DrawText', windowPtr, (txteye), screenCenter(1)/2, screenCenter(2), 0,  backLUT(1));
        Speak(txteye)
        Screen('Flip', windowPtr); % flip to the information screen
        [mx,my,buttons] = GetMouse; % wait for mouse button release before starting experiment
        while any(buttons) % if already down, wait for release
            [mx,my,buttons] = GetMouse(windowPtr);
        end
        while ~any(buttons) % wait for new press
            [mx,my,buttons] = GetMouse(windowPtr);
        end

        nConds=length(blurWidth); % how many letter contrasts to test
        slopeEsts=2*ones(size(blurWidth)); % start with a rough estimate of slope
        minErrEsts=2*ones(size(blurWidth)); % start with a rough estimate of internal error
        startAcuityBnds=[-0.3 1]'; % start with logMAR chart range

        nCells=nRows*nCols;
        [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
        trialSeed=zeros(nTrials, nConds);
        for condNo=1:nConds
            for trialNo=1:nTrials % set up data structure for each chart
                trialRecord(condNo,trialNo) = struct('trialSeed',0,'blurWidth',blurWidth(condNo),'targSizePerLoc',[],'targLocs',[],'targSize', [],'targOri', [],'matchOri', [],'oriErr', [],'respCorrect', [],'stimSeen', [],'chartTime', []); % start empty array
            end
        end
        %% run experiment
        testGapWidthDeg=360/(5*pi); % gap is always 22.9 deg (1/5 of the line width)

        expStart=tic; % start timer
        for trialNo=1:nTrials % work through all trials
            for condNo=1:nConds % work through all trials
                chartStart=tic; % time for each chart
                if trialNo>1 % there has been at least 1 chart, fit all data for next chart
                    testLevel=[]; % set blank list of all levels for this contrast
                    respCorrect=[]; % blank list of all stimuli seen for this contrast
                    for trialSoFar=1:nTrials % work through all trials
                        testLevel=[testLevel trialRecord(trialSoFar).targSizePerLoc]; % concatenate all test sizes
                        %                     angErr=[angErrtrialRecord(condNo,trialNo).oriErr]; % concatenate all orientation errors
                        respCorrect=[respCorrect  trialRecord(trialSoFar).respCorrect]; % concatenate all responses
                    end
                    %                 fitObj=fitPFuncOriErr(tLevel,abs(angErr),1); % fit with current orientation errors
                    fitObj=fitPFuncnAFC_FInDAcuity( testLevel, respCorrect, nAFC, 0); % use Matlab's fit function to fit cumulative Gaussian to the data
                    %                 logThreshEst=fitObj.thresh; % store current threshold estimate
                    acuityBnds = confint(fitObj, 0.95); % 95% confidence intervals on threshold estimate
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

                Screen('Flip', windowPtr); % clear screen

                SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
                trialRecord(condNo,trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
                rng(trialRecord(condNo,trialNo).trialSeed); % seed the random number generator

                trialRecord(condNo,trialNo).targLocs=Shuffle(1:nCells); % pick random target locations for target sizes
                trialRecord(condNo,trialNo).targSize=linspace(minTestLevel, maxTestLevel, nCells); % log spaced range of sizes between upper and lower test sizes in arcmin
                trialRecord(condNo,trialNo).targSizePerLoc(trialRecord(condNo,trialNo).targLocs)=trialRecord(condNo,trialNo).targSize; % fill target locations with corresponding contrast

                trialRecord(condNo,trialNo).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
                trialRecord(condNo,trialNo).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
                trialRecord(condNo,trialNo).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
                trialRecord(condNo,trialNo).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

                xStimCenter=screenCenter(1)+xFactor*(imSize+interCellGapPix); % x center of each cell on threen for this size target
                yStimCenter=screenCenter(2)+yFactor*(imSize+interCellGapPix); % y center of each cell on threen for this size target

                if blurWidth(condNo)>0 % some blur applied
                    myFilter=fspecial('disk',  blurWidth(condNo)*pixPerDeg); % create filter kernel
                end

                for cellNum=1:nCells % draw all random orientation rings
                    gapWidthPix=(10.^trialRecord(condNo,trialNo).targSizePerLoc(cellNum))/pixSizeArcmin; % convert logMAR to arcmin, then to pixels in this system
                    cSize=round(5*gapWidthPix);
                    cImSize=ceil(cSize+4*blurWidth(condNo)*pixPerDeg);
                    mylandoltC=makeLandoltC(cImSize, cSize);
                    mylandoltC(mylandoltC==0)=backLUT;
                    mylandoltC(mylandoltC==1)=targLUT;

                    if blurWidth(condNo)>0 % some blur applied
                        mylandoltC=imfilter(mylandoltC, myFilter, 'replicate');
                    end

                    testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                    trialRecord(condNo,trialNo).targOri(cellNum)=testAngle; % update record for test angle
                    %             destRect=CenterRectOnPoint([0 0 5*gapWidthPix 5*gapWidthPix ], xStimCenter(cellNum), yStimCenter(cellNum)); % rect for target size in this cell
                    destRect=CenterRectOnPoint([0 0 cImSize cImSize], xStimCenter(cellNum), yStimCenter(cellNum)); % rect for target size in this cell

                    myTex=Screen('MakeTexture', windowPtr,mylandoltC);%, [],[],2); % write image to texture, scaled
                    Screen('DrawTexture', windowPtr, myTex, [0 0 cImSize cImSize],destRect,testAngle); % draw texture to on screen location
                    Screen('Close', myTex);

                    %             Screen('FrameArc',windowPtr,targLum,destRect,testAngle+testGapWidthDeg/2,360-testGapWidthDeg,gapWidthPix,gapWidthPix); % arc of required size with required gap
                    destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                    Screen('FrameOval', windowPtr,respLUT,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
                end
                Screen('TextSize',windowPtr,48);                               % put some instructions on screen

                [~, ~, textBox] = DrawFormattedText(windowPtr,  'Click on the orientation of each C:', 100, 100, 0); % give instructions to subject
                % if nTrials~=trialNo
                %       Screen('FillOval',windowPtr,[100 100 100], nextRect);  % draw exit button, grey
                %       DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [128 128 128], [],[],[],[],[],nextRect); % grey 'NEXT' bottom
                %   elseif nTrials==trialNo
                %       Screen('FillOval',windowPtr,[100 100 100], nextRect);  % draw exit button, grey
                %       DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [128 128 128], [],[],[],[],[],nextRect); % grey 'NEXT' bottom
                % end

                Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer

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
                    %             if mx < 1.5*imSize && my < 1.5*imSize % observer clicked finished area
                    if mx > nextRect(1) % observer clicked finished word
                        clickedExit=1;
                    else
                        mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                        [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                        %                 mouseOriWRTChosenBox=90+atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                        mouseOriWRTChosenBox=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                        destRect=CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback
                        Screen('FrameOval', windowPtr,backLUT*0.9,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
                        %                 Screen('FrameArc',windowPtr,[255 200 200],destRect,mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                        Screen('FrameArc',windowPtr,[10 10 10],destRect,90+mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                        Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                    end
                    trialRecord(condNo,trialNo).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
                    trialRecord(condNo,trialNo).oriErr(respNum)=diff(unwrap([trialRecord(condNo,trialNo).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % calculate difference betwen actual and reporte orientation
                    trialRecord(condNo,trialNo).stimSeen(respNum)=1; % mark this cell as response made
                    %computing 'next' buttom
                    if sum( trialRecord(condNo,trialNo).stimSeen)==nCells && nTrials~=trialNo %change color of next button to green once all cells have been clicked, say NEXT
                        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                        DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [128 128 128], [],[],[],[],[],nextRect); %  char(26)
                    elseif  sum( trialRecord(condNo,trialNo).stimSeen)==nCells && nTrials==trialNo %change color of next button to green once all cells have been clicked, say END
                        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button, green
                        DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [128 128 128], [],[],[],[],[],nextRect); %  char(26)
                    end
                    Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer

                    if abs(trialRecord(condNo,trialNo).oriErr(respNum))<(360/(nAFC*2)) % if within nAFC error
                        trialRecord(condNo,trialNo).respCorrect(respNum)=1; % score correct
                    else
                        trialRecord(condNo,trialNo).respCorrect(respNum)=0; % score incorrect
                    end
                end
                trialRecord(condNo,trialNo).chartTime=toc(chartStart);
                if saveImage==0 % saving demo image?
                    myImfileName=sprintf('AIMAcuity%d.jpg', trialNo); % create filename
                    myImage=Screen('GetImage', windowPtr); % grab screnshot
                    imwrite(myImage,myImfileName); % write screenshot to image file
                end
            end % end conds loop
        end % end trials loop
        expDuration=toc(expStart); % stop timer - how long was experiment?
        %Screen('CloseAll'); % close all windows

        save(dataFile,'ExpDate','trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR'); % save parameters

        % threshAcuity=zeros(1, nConds); % set up empty results variable
        % logMAR2Snellen=zeros(1, nConds); % set up empty results variable
        % for condNo=1:nConds
        %     subplot(ceil(nConds^0.5),ceil(nConds/ceil(nConds^0.5)),condNo);
        %     testLevel=[]; % set blank list of all levels for this SF
        %     respCorrect=[]; % blank list of all stimuli seen for this SF
        %     for trialSoFar=1:nTrials % work through all trials
        %         testLevel=[testLevel trialRecord(condNo,trialSoFar).targSizePerLoc]; % concatenate all test sizes
        %         respCorrect=[respCorrect  trialRecord(condNo,trialSoFar).respCorrect]; % concatenate all orientation errors
        %     end
        %     [acuityFitobject,~,threshAcuity(condNo)]=fitPFuncnAFC_FInDAcuity(testLevel, respCorrect, nAFC, 1); % use matlab's fit function to fit cumulative Gaussian to the data
        %     logMAR2Snellen(condNo)=20*(10^threshAcuity(condNo)); % logMAR to Snellen
        %     title([testEye,sprintf(': Participant %s Acuity: logMAR %1.2f Snellen 20/%.1f', sName, threshAcuity(condNo), logMAR2Snellen(condNo))]);
        %     xlabel('Letter Size (logMAR)'); % label the x axis
        %     ylabel('Proportion Correct'); % label y axis
        %     ylim([0 1]); % fix upper and lower bounds
        % end
        %
        % save(dataFile,'logMAR2Snellen','threshAcuity', 'trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect', 'threshAcuity','minLogMAR', 'maxLogMAR');
        %figure();
        angThreshAcuity=zeros(1, nConds);
        anglogMAR2Snellen=zeros(1, nConds);
        for condNo=1:nConds
            subplot(ceil(nConds^0.5),ceil(nConds/ceil(nConds^0.5)),condNo);
            testLevel=[]; % set blank list of all levels for this SF
            respErr=[]; % blank list of all stimuli seen for this SF
            stimSeen=[]; % blank list of all stimuli seen for this SF
            for trialSoFar=1:nTrials % work through all trials
                testLevel=[testLevel trialRecord(condNo,trialSoFar).targSizePerLoc]; % concatenate all test sizes
                respErr=[respErr trialRecord(condNo,trialSoFar).oriErr]; % concatenate all orientation errors
                stimSeen=[stimSeen trialRecord(condNo,trialSoFar).stimSeen]; % concatenate all responses made
            end
            respErr(stimSeen==0)=90; % assign error on ignored cells to 90 deg
            [angleFitobject{condNo},~,angThreshAcuity(condNo),~]=fitPFuncOriErr(testLevel',abs(respErr)',90,0); % fit cumulative Gaussian to the ori error data
            % anglogMAR2Snellen(condNo)=20*(10^angThreshAcuity(condNo)); % logMAR to Snellen
            % title([testEye,sprintf(': Participant %s Acuity: logMAR %1.2f Snellen 20/%.1f', sName, angThreshAcuity(condNo), anglogMAR2Snellen(condNo))]);
            % xlabel('Letter Size (logMAR)'); % label the x axis
            % ylabel('Proportion Correct'); % label y axis
            % ylim([0 180]); % fix upper and lower bounds

            angErr_mean(condNo)= mean(respErr);%mean error (bias)
            angErr_std(condNo)= std(respErr);%std error (bias)

        end
        save(dataFile,'ExpDate','angErr_mean','angErr_std','angleFitobject','angThreshAcuity','trialRecord','expDuration',...
            'sName','correction','blurWidth','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','testEye');
        Screen('Flip', windowPtr); % clear screen
        Screen('DrawText', windowPtr, ('Next eye loading...'), screenCenter(1), screenCenter(2), 0,  backLUT(1));
        Screen('Flip', windowPtr); % clear screen
        close all
    end%ending for randomeye loop
    Screen('CloseAll');

catch exception
    Screen('CloseAll');
    expDuration=toc(expStart); % stop timer - how long was experiment?

    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIM_LowLumVA_BlurFailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sAIM_LowLumVA_BlurfailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'ExpDate','exception','trialRecord','expDuration','tContrast','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect','minLogMAR', 'maxLogMAR');
end




% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/10/2022
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject ] = fitPFuncnAFC_FInDAcuity( level, response, nAFC, graphData )
% Fit data from nAFC using cumulative gaussian psychometric function.
% Function returns fitted parameters, goodness of fit, and 95%
% confidence intervals for the fitted parameters.
% % example data
% level=[0.500000000000000;0.445625469066873;0.397164117362141;0.353972892192069;0.315478672240097;0.281170662595175;0.250593616813636;0.223341796075482;0.199053585276749;0.177406694616788;0.158113883008419;0.140919146563223;0.125594321575479;0.111936056928417;0.0997631157484441;0.125594321575479;0.111936056928417;0.0997631157484441;0.0889139705019463;0.0792446596230558;0.0997631157484441;0.0889139705019463;0.0792446596230558;0.0706268772311378;0.0629462705897085;0.0561009227150983;0.0500000000000001;0.0445625469066874;0.0397164117362142;0.0353972892192070;0.0445625469066874;0.0397164117362142;0.0353972892192070;0.0315478672240097;0.0397164117362142;0.0500000000000001;0.0629462705897085;0.0561009227150983;0.0706268772311379;0.0629462705897085;0.0792446596230559;0.0997631157484443;0.0889139705019464;0.0792446596230559;0.0706268772311380;0.0629462705897086;0.0561009227150984;0.0706268772311380;0.0889139705019465;0.111936056928417];
% response=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;0;1;1;1;1;1;1;1;1;1;0;1;1;1;0;0;0;1;0;1;0;0;1;1;1;1;1;0;0;0;0];
% nAFC=2;
% graphData=1;
if nAFC>0; pGuess=1/nAFC; % guess rate
else pGuess=0;
end
arrangedData=arrangeData(double(level), response);
tLevel=arrangedData(:,1);
nTests=arrangedData(:,2);
pCorrect=arrangedData(:,4);
dataErrEst=arrangedData(:,8);
%   nAFC cumulative gaussian
if any(isnan(dataErrEst))
    dataErrEst=ones(size(dataErrEst));
end
% % ft = fittype( @(thresh, slope,x)(pGuess + (1-pGuess).*(0.5+0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
ft = fittype( @(thresh, slope,x)(pGuess + (1-pGuess).*(0.5+0.5.*erf((x-thresh)./(sqrt(2).*slope)))));
%[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst);
%[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)]);%, 'Weights', 1./dataErrEst);
[fitobject, gof] = fit(tLevel, pCorrect,ft, 'StartPoint', [mean(tLevel) std(tLevel)], 'Weights', 1./dataErrEst,  'Upper',[max(tLevel) range(tLevel)],  'Lower',[min(tLevel) min(diff(tLevel))]);
% [fitobject, gof] = fit(tLevel, pCorrect,ft);%, 'StartPoint', [mean(tLevel) std(tLevel)]);
thresholdEst=fitobject.thresh;
ci = confint(fitobject);
expVals=fitobject(tLevel); % evaluate fit at test levels
chiSqFitReject=chi2gof(sum(((expVals-pCorrect).^2)./expVals), 'Alpha', 0.05, 'NBins',length(tLevel), 'NParams', 2 ); % test hypothesis that fit is significantly different from data (0 means fit is accepted)
if graphData==1
    % ci95 = predint(fitobject,tLevel,0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    %          scatter(tLevel, pCorrect,nTests+20, 'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    %         hold on
    errorbar(tLevel, pCorrect, dataErrEst, dataErrEst, 'bo');%, 'MarkerSize',nTests+2);
    hold on
    ylim([0,1]);
    plot(fitobject, 'r-');
    if ~any(isnan(ci(:))) % no NaNs in confidence intervals
        plot(tLevel,ci,'g--'); % plot 95% confidence interval on fit line
    end
    hold off
    legend HIDE
    box off
    %         set(gca, 'XScale', 'log');
    xlabel('Test Level');
    ylabel('Proportion Correct');
    title(sprintf('Threshold=%.2f, 95%% CI (%.2f,%.2f)', thresholdEst, ci(1,1), ci(2,1)));
end
end

% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/10/2022
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject, lettAcuityEquiv] = fitPFuncOriErr( tSize, mErr,ranErr, graphData,eye,sName, criterion )
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
% [fitobject, gof] = fit(tOri, mErr,ft, 'StartPoint', [mean(tOri) std(tOri) 10], 'Weights', 1./dataErrEst,  'Upper',[max(tOri) range(tOri)],  'Lower',[min(tOri) min(diff(tOri))]);
[fitobject, gof] = fit(tSize, mErr,ft, 'StartPoint', [mean(tSize) std(tSize) 10],  'Upper',[max(tSize) range(tSize) ranErr/2],  'Lower',[min(tSize) 0.05 0]); %Startpoint mean std and arbitary estimate of 10 degree err; max= can't get higher than max size,
% [fitobject, gof] = fit(tLevel, pCorrect,ft);%, 'StartPoint', [mean(tLevel) std(tLevel)]);
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
    if ~any(isnan(ci(:))) % no NaNs in confidence intervals
        ci95 = predint(fitobject,tOriSort, 0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    else
        ci95=nan;
    end
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
    %set(gca, 'XScale', 'log');
    title([sprintf('VA=%.2f (%.2f,%.2f)', thresholdEst, ci(1,1), ci(2,1))],'FontSize',18);
    %         title([sprintf('  %s| ID:%s | Threshold=  %1.2f ', eye,sName,thresholdEst)],'FontSize',18);
    subtitle(['Noise= ' num2str(fitobject.minErr,2), ' \circ  | Slope= ',num2str(fitSlope,2)]); %
    xlabel('Visual Acuity [logMAR]','FontSize',14);
    ylabel('Angular Error [\circ]','FontSize',14);
end
end


function myLandoltC = makeLandoltC(imSize, cSize)
% Make image of Standard Landolt C stimulus
% linewidth is 1/5 optotype size
if nargin==1
    cSize=imSize;
end
[X,Y]=meshgrid(linspace(-imSize/2, imSize/2, imSize),linspace(-imSize/2, imSize/2, imSize)); % 2D matrix of radial distances from centre
radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
myLandoltC=zeros(imSize); % set up background
myLandoltC(radDist<=cSize/2)=1; % make disk of 1's of required size
myLandoltC(radDist<=cSize/2-cSize/5)=0; % convert to annulus
myLandoltC(abs(Y)<imSize/(10*imSize/cSize) & X>0)=0; % % add horizontal right gap
% myLandoltC(abs(X)<imSize/10 & Y<0)=0; % % add vertical top gap
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


