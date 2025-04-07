% AIM BCR Binocular Contrast Ratio
% Program to measure contrast ratio at which interocular contrast is balanced
% uses Ding and Sperling cyclopean phase methiod: Proc Natl Acad Sci U S A. 2006 Jan 24;103(4):1141-6. doi: 10.1073/pnas.0509629103. Epub 2006 Jan 12.
% Presents sequence of grids containing random orientation blurred edge stimuli of a range of interocular contrasts
% Contrast ratio based on threshold estimate -2 std dev to + 2 stdevs
% Observer clicks perceived orientation of each ring
% Responses proportional to contrast gain in each eye
% fit with cumulative gaussian, mean and stnadard error defines size range for next chart

% 2024  PJB
% WhiteVal=252.5; % luminance of LUT 255
% BlackVal=1.5; % luminance of LUT 0
% meanLUT=LumDesired/WhiteVal*255;

% This software is patented and owned by Northeastern University, Boston,
% USA; exclusively liscenced by PerZeption Inc.
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021

%% Confidential and exclusively use by Prof Zahide for research purposes only. Not for third parties nor commercial applications; not for diagnosis nor screening.

clear all; close all; commandwindow;n=1; % close windows and type into command window

whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% read in testing parameters:
prompt = {'Subject Initials', 'Stereo Method','Background level (0-1)', 'Ring Diameter (deg)','Blur Stdev (list, deg)','Cell Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)','With(1) of Without(0) Correction'};
dlg_title = 'AIM BCR';
num_lines = 1;
def = {'XX', '-1', '0.5', '6','8', '0.5', '4','4','3', '120', '60','1'}; % some defaul values
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % assign experimenter values to experimnet
stereoMode=str2num(char(answer(2,1))); % method of dichoptic presentation
meanLUM=str2num(char(answer(3,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
targSizeDeg=str2num(char(answer(4,1))); % Size of C - Pelli Robson letters are 4.8cm, @1m = 2.75 deg
cellSizeDeg=str2num(char(answer(5,1))); % diameterof cell rings
blurStdDeg=str2num(char(answer(6,1))); % std of blur in degree
nRows=str2num(char(answer(7,1))); % $ rows and columns on each chart
nCols=str2num(char(answer(8,1)));
nTrials=str2num(char(answer(9,1))); % # charts to run
scrnWidthCm=str2num(char(answer(10,1))); % screen width (cm)
viewDistCm=str2num(char(answer(11,1))); % observer's distance from screen
correction=str2num(char(answer(12,1)));

rng('default'); % seed random number generator
gammaVal=2.0; % gamma on this system
LMean=255*meanLUM.^(1/gammaVal); % Gamma-corrected background

nAFC=0; % angle for error scoring
saveImage=0; % save stimuli or not
interCellGapDeg=2;

try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    PsychImaging('PrepareConfiguration');
    if stereoMode<0
        PsychImaging('AddTask', 'General', 'InterleavedLineStereo', 0);
        [windowPtr, winRect] = PsychImaging('OpenWindow', whichScreen, LMean);    % open experimental window with required parameters
        aspectRatio=0.5; % correct for doubling image size
    else
        [windowPtr, winRect] = PsychImaging('OpenWindow', whichScreen, LMean, [],  [], [], stereoMode);%, [0 0 3840 2160]); % open the window for the experiment
        aspectRatio=1.0; % correct for doubling image size
    end
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    cellDiamPix=round(cellSizeDeg*pixPerDeg); % size of each cell in pixels on this system
    responseRingPenWidth=pixPerDeg*0.5; % line width for rresponse ring - TODO justify 0.25 deg
    interCellGapPix=ceil(interCellGapDeg*pixPerDeg);
    targSizePix=round(targSizeDeg*pixPerDeg); % size of target on this system
    testGapWidthDeg=360/(5*pi); % gap is always 22.9 deg (1/5 of the line width)
    minCRange=0.2; % minium range o contrasts on screen
    nextSizeDeg=cellSizeDeg/2; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg aspectRatio*nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left

    blurStdevPix=blurStdDeg*pixPerDeg;
    % oriOffset=360*(blurStdevPix/(2*pi()*targSizePix/2)); % # degrees for 1 standard deviation at half radius

    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,32);                               % put some instructions on screen
    textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz \\nClick on the orientation of each ring \\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate);

    Screen('SelectStereoDrawBuffer', windowPtr, 0); % draw left eye stimulus
    DrawFormattedText(windowPtr, textToObserver, 100, 100, 0, LMean); % put this message on screen
    Screen('SelectStereoDrawBuffer', windowPtr, 1); % draw right eye stimulus
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

    nBlurs=length(blurStdDeg); % how many letter contrasts to test
    slopeEsts=2*ones(size(blurStdDeg)); % start with a rough estimate of slope
    minErrEsts=2*ones(size(blurStdDeg)); % start with a rough estimate of internal error
    startContrastBnds=[0 1]'; % OS 0% to OS 100%

    [X,Y]=meshgrid(linspace(-cellDiamPix/2, cellDiamPix/2, cellDiamPix),linspace(-cellDiamPix/2, cellDiamPix/2, cellDiamPix)); % 2D matrix of radial distances from centre
    radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
    respAnnulusIm=zeros(cellDiamPix); % set up background
    respAnnulusIm(radDist<=cellDiamPix/2)=1; % make disk of 1's of required size
    respAnnulusIm(radDist<=cellDiamPix/2-responseRingPenWidth)=0; % convert to annulus
    angIm=atan2(Y, -X)/pi;  % orientation of each pixel, scaled -1 to +1
    responseRing=ones([cellDiamPix cellDiamPix 2]); % set up 2 layer image for window and alpha
    responseRing(:,:,1)=255*(0.5+0.5*angIm).^(1/gammaVal); % gamma corrected ramp
    responseRing(:,:,2)=255*(respAnnulusIm);
    responseRingRect=[0 0 cellDiamPix cellDiamPix]; % dimensions of response ring source
    responseRingOnScreenRect=[0 0 cellDiamPix cellDiamPix*aspectRatio]; % dimensions of response ring on screen

    srcRect=[0 0 targSizePix targSizePix]; % rect for this sized cell
    onScreenRect=[0 0 targSizePix targSizePix*aspectRatio]; % onscreen rect for targets

    nCells=nRows*nCols;
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
    trialSeed=zeros(nTrials, nBlurs);

    for trialNo=1:nTrials % set up data structure for each chart
        for blurNum=1:nBlurs % and for each contrast
            trialRecord(trialNo,blurNum) = struct('trialSeed',0,'targContrastPerLoc',[],'targLocs',[],'targContrast', [],'targOri', [],'matchOri', [],'oriErr', [],'stimSeen', [],'chartTime', []); % start empty array
        end
    end
    % if stereoMode==6 || correction==1
    %     Speak('Please wear the red blue glasses with your glasses or contact lenses.')
    % elseif stereoMode==6 || correction==0
    %     Speak('Please wear the red blue glasses.')
    % else
    % end
    %intrcutions
            Speak('Click where the light and dark half of each ring are starting from.')

    %% run experiment

    expStart=tic; % start timer
    for trialNo=1:nTrials % work through all trials

        for blurNum=1:nBlurs % work though the list of test contrasts
            mouseX=[]; mouseY=[]; mouseclick=[];%empty variables for mouse data, JSk edited 2024

            chartStart=tic; % time for each chart
            if trialNo>1 % there has been at least 1 chart, fit all data for next chart
                angErr=[]; % set blank list of all levels for this contrast
                testLevel=[]; % set blank list of all levels for this contrast
                stimSeen=[]; % blank list of all stimuli seen for this SF
                for trialSoFar=1:nTrials % work through all trials
                    testLevel=[testLevel trialRecord(trialSoFar,blurNum).targContrastPerLoc]; % concatenate all test sizes
                    angErr=[angErr trialRecord(trialSoFar,blurNum).oriErr]; % concatenate all orientation errors
                    stimSeen=[stimSeen trialRecord(trialSoFar,blurNum).stimSeen]; % concatenate all responses made
                end
                seenTestLevel=testLevel(stimSeen==1); % only process stimuli with response -level
                seenAngErr=angErr(stimSeen==1); % only process stimuli with response - report error
                fitObj=fit_BCR_OriErrorFunc( seenTestLevel', seenAngErr', 0 , -1);
                csBnds = confint(fitObj, 0.95); % 95% confidence intervals on contrast threshold estimate
                if any(isnan(csBnds)) % fit failed to estimate acuity and 95% CIs
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

            baseRing=makeBCR_Ring(targSizePix, blurStdevPix(blurNum));
            oriOffset=360*(blurStdevPix(blurNum)/(2*pi()*targSizePix/2)); % # degrees for 1 standard deviation at half radius

            minTestLevel=max([0 contrastBound(1)]); % lower 99% CI contrast, minimum 0%
            maxTestLevel=min([1 contrastBound(2)]); % uppper 99% CI contrast, maximum 100%
            if (maxTestLevel-minTestLevel) < minCRange % less than 2dB range
                midTestLevel=mean([minTestLevel maxTestLevel]);
                minTestLevel=midTestLevel-minCRange/2; % lower min to span min range
                maxTestLevel=midTestLevel+minCRange/2; % raise max to span min range
            end
            Screen('Flip', windowPtr); % clear screen

            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(trialNo, blurNum).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(trialNo, blurNum).trialSeed); % seed the random number generator

            trialRecord(trialNo,blurNum).targLocs=Shuffle(1:nCells); % pick random target locations for target contrasts
            trialRecord(trialNo,blurNum).targContrast=linspace(minTestLevel, maxTestLevel, nCells); % log spaced range of contrasts between upper and lower test contrasts
            trialRecord(trialNo,blurNum).targContrastPerLoc(trialRecord(trialNo,blurNum).targLocs)=trialRecord(trialNo,blurNum).targContrast; % fill target locations with corresponding contrast

            trialRecord(trialNo,blurNum).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,blurNum).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
            trialRecord(trialNo,blurNum).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect

            xStimCenter=screenCenter(1)+xFactor*(cellDiamPix+interCellGapPix); % x center of each cell on threen for this size target
            yStimCenter=screenCenter(2)+yFactor*aspectRatio*(cellDiamPix+interCellGapPix); % y center of each cell on threen for this size target

            for cellNum=1:nCells % draw all random orientation rings
                testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
                trialRecord(trialNo,blurNum).targOri(cellNum)=testAngle; % update record for test angle
                destRect=CenterRectOnPoint(onScreenRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring

                myTarg=meanLUM+meanLUM*baseRing*trialRecord(trialNo,blurNum).targContrastPerLoc(cellNum); % scale to mean luminance - contrast * C
                myTarg=255*myTarg.^(1/gammaVal); % gamma correct
                myTarg=imrotate(myTarg,testAngle-oriOffset,'nearest', 'crop'); % left eye CCW rorated 1 std
                myTarg(myTarg==0)=LMean;
                myTex=Screen('MakeTexture', windowPtr,myTarg);%, [],[],2); % write image to texture, scaled
                Screen('SelectStereoDrawBuffer', windowPtr, 0); % draw left eye stimulus
                % Screen('DrawTexture', windowPtr, myTex, srcRect,destRect,testAngle-oriOffset); % draw texture to on screen location
                Screen('DrawTexture', windowPtr, myTex, srcRect,destRect); % draw texture to on screen location
                Screen('Close', myTex);

                myTarg=meanLUM+meanLUM*baseRing*(1-trialRecord(trialNo,blurNum).targContrastPerLoc(cellNum)); % scale to mean luminance - contrast * C
                myTarg=255*myTarg.^(1/gammaVal); % gamma correct
                myTarg=imrotate(myTarg,testAngle+oriOffset,'nearest', 'crop'); % left eye CCW rorated 1 std
                myTarg(myTarg==0)=LMean;
                myTex=Screen('MakeTexture', windowPtr,myTarg);%, [],[],2); % write image to texture, scaled
                Screen('SelectStereoDrawBuffer', windowPtr, 1); % draw right eye stimulus
                Screen('DrawTexture', windowPtr, myTex, srcRect,destRect); % draw texture to on screen location
                Screen('Close', myTex);

                destRect=CenterRectOnPoint(responseRingOnScreenRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
                respRingOri=360*rand();
                myTarg=imrotate(responseRing,respRingOri,'nearest', 'crop'); % left eye CCW rorated 1 std
                myTex=Screen('MakeTexture', windowPtr, myTarg); % make alpha image as texture
                Screen('SelectStereoDrawBuffer', windowPtr, 0); % draw left eye stimulus
                %Screen('DrawTexture', windowPtr, myTex, responseRingRect, destRect); % draw OS response ring
                Screen('FrameOval', windowPtr,LMean ,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells

                Screen('SelectStereoDrawBuffer', windowPtr, 1); % draw right eye stimulus
                %Screen('DrawTexture', windowPtr, myTex, responseRingRect, destRect); % draw OD response ring
                Screen('FrameOval', windowPtr,LMean ,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells

                Screen('Close', myTex);
            end

            % Screen('SelectStereoDrawBuffer', windowPtr, 0); % draw left eye stimulus
            % Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
            % DrawFormattedText(windowPtr, sprintf('%s', 'next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
            % Screen('SelectStereoDrawBuffer', windowPtr, 1); % draw right eye stimulus
            % Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
            % DrawFormattedText(windowPtr, sprintf('%s', 'next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)

            Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer

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
                %                 if mx < textBox(3) && my < textBox(4) % observer clicked finished area
                if mx > nextRect(1) % observer clicked finished word
                    clickedExit=1;
                else
                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my*aspectRatio-yStimCenter).^2); % calculate distance from each box
                    [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                    mouseOriWRTChosenBox=atan2d(-(my*aspectRatio-yStimCenter(respNum))/aspectRatio,mx-xStimCenter(respNum)); % perceived orientation for this response
                    destRect=CenterRectOnPoint(responseRingOnScreenRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback

                    respRingOri=mouseOriWRTChosenBox;
                    myTarg=imrotate(responseRing,respRingOri,'nearest', 'crop'); % left eye CCW rorated 1 std
                    myTex=Screen('MakeTexture', windowPtr, myTarg); % make alpha image as texture
                    Screen('SelectStereoDrawBuffer', windowPtr, 0); % draw left eye stimulus
                    Screen('DrawTexture', windowPtr, myTex, responseRingRect, destRect); % draw OS response ring
                    Screen('SelectStereoDrawBuffer', windowPtr, 1); % draw right eye stimulus
                    Screen('DrawTexture', windowPtr, myTex, responseRingRect, destRect); % draw OD response ring
                    Screen('Close', myTex);

                    Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
                end
                trialRecord(trialNo,blurNum).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
                trialRecord(trialNo,blurNum).oriErr(respNum)=diff(unwrap([trialRecord(trialNo,blurNum).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % calculate difference betwen actual and reported orientation
                trialRecord(trialNo,blurNum).stimSeen(respNum)=1; % mark this cell as response made
                %edited by JSk,2024
                %collect mouse data
                trialRecord( trialNo,blurNum).MouseX=mouseX; %x coordinates for mouse
                trialRecord( trialNo,blurNum).MouseY=mouseY; %y coordinates for mouse
                trialRecord( trialNo,blurNum).Mousebuttons=mouseclick; %x coordinates for mouse
                %computing 'next' buttom
                if sum( trialRecord(trialNo,blurNum).stimSeen)==nCells && nTrials~=trialNo %change color of next button to green once all cells have been clicked, say NEXT
                    Screen('FillOval',windowPtr,[100 100 100], nextRect);  % draw exit button, green
                    DrawFormattedText(windowPtr, sprintf('%s', 'NEXT'), 'center', 'center', [128 128 128], [],[],[],[],[],nextRect); %  char(26)
                elseif  sum( trialRecord(trialNo,blurNum).stimSeen)==nCells && nTrials==trialNo %change color of next button to green once all cells have been clicked, say END
                    Screen('FillOval',windowPtr,[100 100 100], nextRect);  % draw exit button, green
                    DrawFormattedText(windowPtr, sprintf('%s', 'END'), 'center', 'center', [128 128 128], [],[],[],[],[],nextRect); %  char(26)
                end
                Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
            end
            trialRecord(trialNo,blurNum).chartTime=toc(chartStart);
        end % end SFs loop
        if saveImage ==1% saving demo image?
            myImfileName=sprintf('AIM_BCR%d.jpg', trialNo); % create filename
            myImage=Screen('GetImage', windowPtr); % grab screnshot
            imwrite(myImage,myImfileName); % write screenshot to image file
        end
    end % end trials loop
    expDuration=toc(expStart); % stop timer - how long was experiment?
    Screen('CloseAll'); % close all windows
    %% save data
    if correction ==1
        correction='wCorr';
    elseif correction ==0
        correction='woCorr';
    end
    testSName=[sName, '_' correction '_' num2str(n)];

    while exist([testSName,'AIM_BCR.mat'],'file') ~= 0 % check if data file exists with this subject ID
        n=n+1;
        testSName=[sName, '_' correction '_' num2str(n)];
    end
    dataFile=sprintf('%sAIM_BCR.mat', testSName);  % file path to unique Mat files
    save(dataFile,'trialRecord', 'expDuration','sName', 'stereoMode', 'meanLUM','targSizeDeg','cellSizeDeg','blurStdDeg', 'nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

   

    balPoint = zeros(1, nBlurs); % set up empty results variable
    for blurNum = 1:nBlurs % plot figure with fits for interleaved conditions - attempt to fit all data
        %subplot(ceil(sqrt(nBlurs)), ceil(nBlurs/ceil(sqrt(nBlurs))), blurNum);
        figure
        angErr = []; % set blank list of all levels for this contrast
        testLevel = []; % set blank list of all levels for this contrast
        stimSeen = []; % blank list of all stimuli seen for this SF
        for trialSoFar = 1:nTrials % work through all trials
            testLevel = [testLevel, trialRecord(trialSoFar, blurNum).targContrastPerLoc]; % concatenate all test sizes
            angErr = [angErr, trialRecord(trialSoFar, blurNum).oriErr]; % concatenate all orientation errors
            stimSeen = [stimSeen, trialRecord(trialSoFar, blurNum).stimSeen]; % concatenate all responses made
        end
        seenTestLevel = testLevel(stimSeen == 1); % only process stimuli with response - level
        seenAngErr = angErr(stimSeen == 1); % only process stimuli with response - report error
        [fitObj, ~, balPoint(blurNum)] = fit_BCR_OriErrorFunc(seenTestLevel', seenAngErr', 1, -1);
        subtitle(['Condition: ' num2str(blurStdDeg(blurNum)) ' \circ blur SD'])
        ylim([-185 185])
        line([0 1], [0 0], 'Color', 'black', 'LineStyle', '--'); % horizontal orientation
        line([0.5 0.5], [-185 185], 'Color', 'black', 'LineStyle', '--'); % vertical orientation
        text(0.7, 60, ['Slope: ' num2str(fitObj.slope, 2)])
        text(0.7, 80, ['Scale: ' num2str(fitObj.scale, 2)])

        % Create the inset bar plot
        if balPoint(blurNum)==0.5 % equilibrium
            right_eye_data=balPoint(blurNum)*100;
            left_eye_data=balPoint(blurNum)*100;
        elseif balPoint(blurNum)>0.5 % right eye dominance-blue right-red left
            right_eye_data=abs(1-balPoint(blurNum))*100;
            left_eye_data=balPoint(blurNum)*100;
        elseif balPoint(blurNum)<0.5 % left eye dominance
            right_eye_data=balPoint(blurNum)*100;
            left_eye_data=abs(1-balPoint(blurNum))*100;
        end
        % Sample data for right and left eye
      
        inset_ax = axes('Position', [0.75, 0.25, 0.2, 0.2],'box','off'); % Adjust position as needed
        bar(inset_ax, [1, 2], [(right_eye_data), (left_eye_data)], 'FaceColor', 'flat');
        inset_ax.YLim = [0 100];
        inset_ax.XTickLabel = {'OD', 'OS'};
        ylabel(inset_ax, '[%]');
    end


    save(dataFile,'trialRecord','fitObj', 'balPoint', 'expDuration','sName', 'stereoMode', 'meanLUM','targSizeDeg','cellSizeDeg','blurStdDeg', 'nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters

catch exception
    Screen('CloseAll');
    testSName=[sName, '_' correction '_' num2str(n)];
    while exist([testSName,'AIM_BCR_FailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName, '_' correction '_' num2str(n)];
    end
    dataFile=sprintf('%sAIM_BCR_FailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'exception','trialRecord','expDuration','sName', 'stereoMode', 'meanLUM','targSizeDeg','cellSizeDeg','blurStdDeg', 'nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect'); % save parameters
end



function [ fitobject, gof, thresholdEst, ci ] = fit_BCR_OriErrorFunc(tLevel, oriErr, graphData, slopeSign )

% Fit data from contrast ratio task using cumulative gaussian psychometric function.
% Function returns fitted parameters, goodness of fit, and 95%
% confidence intervals for the fitted parameters.
% % example data
% thresh=0.75;
% slope=0.1;
% scale=16;
% tLevel=linspace(0,1,20);
% oriErr=-(scale.*erf((tLevel-thresh)./(sqrt(2).*slope)))+(scale/8)*randn(size(tLevel));
% fit_BCR_OriErrorFunc( tLevel', oriErr', 1 );

if nargin==3
    slopeSign=1;
end
ft = fittype( @(thresh, slope,scale, x)(scale.*erf((x-thresh)./(sqrt(2).*slope*slopeSign))));

[fitobject, gof] = fit(tLevel, oriErr,ft, 'StartPoint', [mean(tLevel) std(tLevel) max(abs(oriErr))], 'Upper',[max(tLevel) range(oriErr)/2 2*max(abs(oriErr))],  'Lower',[min(tLevel) min(diff(tLevel)) 0.1*max(abs(oriErr))]);

thresholdEst=fitobject.thresh;
ci = confint(fitobject);

if graphData
    xEval=linspace(min(tLevel), max(tLevel));
    ci95 = predint(fitobject,xEval,0.95,'functional','off'); % pull 95% confdence intervals from fitParams at each test level
    scatter(tLevel, oriErr, 'bo');
    hold on
    plot(fitobject, 'r-');
    if ~any(isnan(ci95(:))) % no NaNs in confidence intervals
        plot(xEval,ci95,'g--'); % plot 95% confidence interval on fit line
    end
    hold off
    legend HIDE
    box off
    xlabel('Contrast Ratio');
    ylabel('Report Error [\circ]');
    title(sprintf('Balance point=%.2f, 95%% CI (%.2f,%.2f)', thresholdEst, ci(1,1), ci(2,1)));

end
end


%
% % function myBCRRing = makeBCR_ResponseRing(imSize, blurStdev)
% % % Make image of Standard Landolt C stimulus
% % % linewidth is 1/5 optotype size
%
% imSize=256;
% blurStdev=32;
%
% [X,Y]=meshgrid(linspace(-imSize/2, imSize/2, imSize),linspace(-imSize/2, imSize/2, imSize)); % 2D matrix of radial distances from centre
% radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
%
%
% respAnnulusIm=zeros(imSize); % set up background
% respAnnulusIm(radDist<=imSize/2)=1; % make disk of 1's of required size
% respAnnulusIm(radDist<=imSize/2-8)=0; % convert to annulus
%
% targAnnulusIm=zeros(imSize); % set up background
% targAnnulusIm(radDist<=imSize/2-32)=1; % make disk of 1's of required size
% targAnnulusIm(radDist<=imSize/2-(32+48))=0; % convert to annulus
%
%
% angIm=atan2(Y, -X)/pi;  % orientation of each pixel, scaled -1 to +1
% myBCRRing=(respAnnulusIm.*angIm)+(targAnnulusIm.*imgaussfilt(angIm, blurStdev));
%
% imagesc(myBCRRing)
% colormap gray
% imwrite(uint8(127+127*myBCRRing), 'bcrStim.jpg')
% % end


function myBCRRing = makeBCR_Ring(imSize, blurStdev)
% Make image of Standard Landolt C stimulus
% linewidth is 1/5 optotype size

[X,Y]=meshgrid(linspace(-imSize/2, imSize/2, imSize),linspace(-imSize/2, imSize/2, imSize)); % 2D matrix of radial distances from centre
radDist=(X.^2+Y.^2).^0.5; % radial distance from centre
annulusIm=zeros(imSize); % set up background
annulusIm(radDist<=imSize/2)=1; % make disk of 1's of required size
annulusIm(radDist<=imSize/2-imSize/4)=0; % convert to annulus

angIm=atan2(Y, -X)/pi;  % orientation of each pixel, scaled -1 to +1
myBCRRing=annulusIm.*imgaussfilt(angIm, blurStdev);
end

