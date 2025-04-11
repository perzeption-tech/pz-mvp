% AIM VA
% Program to measure visual acuity based on orientation identification of rotated Landolt C stimulus
% Presents sequence of grids containing random orientation C stimuli of a  range of sizes
% Size of stimuli based on threshold estimate -2 std dev to + 4 stdevs
% Observer clicks perceived orientation of each C and can change response or leave cells empty - scored as incorrect with random orientation error (90 deg)
% Responses classified as 8AFC and fit with cumulative gaussian, mean and stnadard error defines size range for nxt chart
% Data also fit with error function
% 2022  PJB
clear all; close all; commandwindow; % clear all memory, close windows and type into command window
whichScreen=max(Screen('Screens')); % use highest display # for stimuli
ExpDate = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm')
% read in testing parameters:
prompt = {'Subject Initials', 'Background Luminance (0-255)','C Luminance (0-255)','Screen Width (cm)','Viewing Distance (cm)','Numbers of cells'};
dlg_title = 'AIM VA Low Vision';
num_lines = 1;
def = {'XX', '180','0', '70', '100','32'}; % some defaul values
timestamp_experiment_start = datetime('now');
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1)); % assign experimenter values to experimnet
backLUT=str2num(char(answer(2,1))); % Weber contrast of letters
targLUT=str2num(char(answer(3,1))); % Weber contrast of letters
scrnWidthCm=str2num(char(answer(4,1))); % screen width (cm)
viewDistCm=str2num(char(answer(5,1))); % observer's distance from screen
minNumReports=str2num(char(answer(6,1))); % observer's distance from screen

testSName=sName;    % modify sName if subject already exists
n=1;
while exist([testSName,'AIM_VA_LowVision.mat'],'file') ~= 0 % check if data file exists with this subject ID
    n=n+1;
    testSName=[sName,num2str(n)]; % if so create new name until unique name found
end
dataFile=sprintf('%sAIM_VA_LowVision.mat', testSName);  % file path to unique Mat files

rng('default'); % seed random number generator

nAFC=8; % angle for error scoring
saveImage=0; % save stimuli or not
minCellSizeDeg=6; %
gammaVal=2.3; % gamma on this system
% backLUT=255*meanLUM.^(1/gammaVal); % Gamma-corrected background
backLUT =255*(backLUT./255).^(1/gammaVal); % Gamma-corrected background


%respLUT=[backLUT backLUT 255];% response ring color pre-indication
respLUT=[backLUT+10 backLUT+10 backLUT+10];% response ring color pre-indication
respARCLUT=[200 255 200];% response arc color
resppostIndicationLUT=[backLUT-10 backLUT-10 backLUT-10];% response ring color post indication
ranreporterror=90;%random report in C task
try
    Screen('Preference','SkipSyncTests', 1); % PTB hack
    [windowPtr, winRect]=Screen('OpenWindow', whichScreen, backLUT);%[0 0 2000 1800]
    
    frameRate=Screen('FrameRate', windowPtr); % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    scrnHeightDeg=scrnHeightPix/scrnWidthPix*scrnWidthDeg; % height of screen in deg
    
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    nextSize=round(pixPerDeg*1/2); % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSize nextSize],winRect(3)-nextSize/2,screenCenter(2)); % halfway down the screen, away from the left
    %edited JSk, 2024
    %increase gap between cells
    interCellGapPix=pixPerDeg*0.5;   % interCellGapPix=0.5*pixPerDeg;
    
    pixSizeArcmin=60/pixPerDeg; % size of eacdh pixel in arcmin
  
    Screen('TextFont',windowPtr,'Arial');
    Screen('TextSize',windowPtr,48);                               % put some instructions on screen
    textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz\\nClick mouse to start', winRect(3)-winRect(1),winRect(4)-winRect(2),frameRate);
    DrawFormattedText(windowPtr, textToObserver, 100, 100, 0, backLUT); % put this message on screen
    Screen('Flip', windowPtr); % flip to the information screen
    [mx,my,buttons] = GetMouse; % wait for mouse button release before starting experiment
    while any(buttons) % if already down, wait for release
        [mx,my,buttons] = GetMouse(windowPtr);
    end
    while ~any(buttons) % wait for new press
        [mx,my,buttons] = GetMouse(windowPtr);
    end
    nTrials=10; % start with more charts than will be needed
    for trialNo=1:nTrials % set up data structure for each chart
        trialRecord(trialNo) = struct('trialSeed',0,'nRows',0,'nCols',0,'targSizePerLoc',[],'targLocs',[],'targSize', [],'targOri', [],'matchOri', [],'oriErr', [],'respCorrect', [],'stimSeen', [],'chartTime', [],'MouseX',[],'MouseY',[],'Mousebuttons',[],'mylandoltC',[]); % start empty array
    end
    
    %% run experiment
    testGapWidthDeg=10; % response gap is 10 deg
    maxPossTargSize=0.65*scrnHeightDeg; % start with full size letter; largest letter is 50% of screen height in degree, ring is 25% in each direction
    
    minPossTargSize=5/pixPerDeg; % minimum target is 5 pixels
    totalNTrials=0;
    trialNo=1;
    nCells_count=1;%it always starts with one single stimulus ;made for cell counter to detect the first 3*5 grid
    maxLogMAR=log10((pixSizeArcmin*(maxPossTargSize*60))/12); % maximum target size logMAR
    minLogMAR=log10(pixSizeArcmin); % minimum target size logMAR
    expStart=tic; % start timer
    % for trialNo=1:nTrials % work through all trials
    while totalNTrials<minNumReports % complete minimum # reports
        chartStart=tic; % time for each char
        mouseX=[]; mouseY=[]; mouseclick=[];%empty variables for mouse data, JSk edited 2024
        
        if trialNo>1 % compute maxtarg size based on previous trials
            testLevel=[]; % set blank list of all levels for this contrast
            angErr=[]; % blank list of all stimuli seen for this contrast
            respCorrect=[]; % blank list of all stimuli seen for this contrast
            for trialSoFar=1:nTrials % work through all trials
                testLevel=[testLevel trialRecord(trialSoFar).targSizePerLoc]; % concatenate all test sizes
                angErr=[angErr trialRecord(trialSoFar).oriErr]; % concatenate all orientation errors
                respCorrect=[respCorrect  trialRecord(trialSoFar).respCorrect]; % concatenate all responses
            end
            
            totalNTrials=length(testLevel); % total # observations so far
            % if totalNTrials>3 % enought to fit
            %             if sum(respCorrect)<totalNTrials && totalNTrials>3 % some errors and enough parameters to fit
            %                 fitObj=fitPFuncnAFC_AIMLOWVISIONAcuity(testLevel, respCorrect, nAFC, 0); % use Matlab's fit function to fit cumulative Gaussian to the data
            %                 maxTargSize=min([maxPossTargSize fitObj.thresh*2]); % maximum value at 2.0 * threshold estimate
            %                 minTargSize=max([minPossTargSize fitObj.thresh/2]); % maximum value at 2.0 * threshold estimate
            if sum(respCorrect)==totalNTrials && totalNTrials<=3   % all correct and not enough data for a 3 parameter fit
                maxTargSize=min(testLevel)/(10.^0.1); % max decrease 0.1 logMAR
                % minTargSize=min(testLevel)/(10.^0.1).^2; % max decrease 0.2 logMAR
                minTargSize=0; %flag value
            elseif sum(respCorrect)<totalNTrials && totalNTrials<=3   % some errors and not enough data for a 3 parameter fit
                
                maxTargSize=min(testLevel)*(10.^0.1); % max increase 0.1 logMAR
                %       % minTargSize=min(testLevel)/(10.^0.1).^2; % max decrease 0.2 logMAR
                % minTargSize=0; %flag value
                if maxPossTargSize< maxTargSize
                    maxTargSize=maxPossTargSize;
                end
            elseif nCells==1 && totalNTrials>3   %  previously only one stimulus presented, and more than 3 trials
                
                if  trialRecord(trialNo-1).respCorrect ==0 %incorrect
                    maxTargSize=min(testLevel)*(10.^0.1); % max increase 0.1 logMAR
                    %       % minTargSize=min(testLevel)/(10.^0.1).^2; % max decrease 0.2 logMAR
                    % minTargSize=0; %flag value
                    if maxPossTargSize< maxTargSize
                        maxTargSize=maxPossTargSize;
                    end
                elseif trialRecord(trialNo-1).respCorrect ==1 %correct
                    maxTargSize=min(testLevel)/(10.^0.1); % max decrease 0.1 logMAR
                    % minTargSize=min(testLevel)/(10.^0.1).^2; % max decrease 0.2 logMAR
                    minTargSize=0; %flag value
                end
                
            elseif   nCells>1 && totalNTrials>3 % more than one stimulus presented previously, and more than 3 trials
                fitObj=fitPFuncnAFC_AIMLOWVISIONAcuity(testLevel, respCorrect, nAFC, 0); % use Matlab's fit function to fit cumulative Gaussian to the data
                maxTargSize=min([maxPossTargSize fitObj.thresh*2]); % maximum value at 2.0 * threshold estimate
                minTargSize=max([minPossTargSize fitObj.thresh/2]); % maximum value at 2.0 * threshold estimate
                
                %             else
                %                 maxTargSize=min(testLevel)/(10.^0.1); % max decrease 0.1 logMAR
                %                 % minTargSize=min(testLevel)/(10.^0.1).^2; % max decrease 0.2 logMAR
                %                 minTargSize=0;
            end
            
            
        else % only for first chart
            threshSize=maxPossTargSize; % start with maximum that can fit on screen
            maxTargSize=threshSize; % maximum value = threshold estimate
            minTargSize=threshSize; % minimum value = threshold estimate
        end
        %edited by JSk 2024,
        if isempty(floor(scrnHeightDeg/(2*maxTargSize)))==1 % check whether at least one single stimulus can be presented (1 row and 1 column)
            trialRecord(trialNo).nRows=1; % how many rows can be fit on screen, up to 4
            trialRecord(trialNo).nCols=1; % how many rows can be fit on screen, up to 4
        elseif (floor(scrnHeightDeg/(2*maxTargSize)))<=0 % check whether at least one single stimulus can be presented (1 row and 1 column)
            trialRecord(trialNo).nRows=1; % how many rows can be fit on screen, up to 4
            trialRecord(trialNo).nCols=1; % how many rows can be fit on screen, up to 4
        else
            trialRecord(trialNo).nRows=min([3 floor(scrnHeightDeg/(2*maxTargSize))]); % how many rows can be fit on screen, up to 4
            trialRecord(trialNo).nCols=min([5 floor(scrnWidthDeg/(2*maxTargSize))]); % how many rows can be fit on screen, up to 4
        end
        nCells=trialRecord(trialNo).nRows*trialRecord(trialNo).nCols;
        
        
        [xFactor, yFactor]=meshgrid((1:trialRecord(trialNo).nCols)-trialRecord(trialNo).nCols/2-0.5,(1:trialRecord(trialNo).nRows)-trialRecord(trialNo).nRows/2-0.5); % relative distance from center of screen in multiples of stim size
        
        
        SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
        trialRecord(trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
        rng(trialRecord(trialNo).trialSeed); % seed the random number generator
        
        if minTargSize==0 % flag to calculate
            minTargSize=maxTargSize/(10.^0.1).^nCells; % max decrease 0.1 logMAR
            %JSk edit, 2024, January
            if minTargSize<minPossTargSize %ensure that size is never smaller than min of screen resolution
                minTargSize=minPossTargSize;
            end
        end
        trialRecord(trialNo).targLocs=Shuffle(1:nCells); % pick random target locations for target sizes
        if nCells==15 &&  isnan(nCells_count(end))==0%when the first time max number of grid appears show smallest possible stimulus
            trialRecord(trialNo).targSize=logspace(log10(minPossTargSize), log10(maxTargSize), nCells);
            nCells_count=[nCells_count; nan];%deleted variable as it served its purpose
        elseif nCells==15 &&  isnan(nCells_count(end))==1%otherwise show + and - 2 STD around the threshold
            trialRecord(trialNo).targSize=logspace(log10(minTargSize), log10(maxTargSize), nCells); % log spaced range of sizes between upper and lower test sizes in arcmin
                   nCells_count=[nCells_count ;0];%
        else % show + and - 2 STD around the threshold, or single stimuli
            trialRecord(trialNo).targSize=logspace(log10(minTargSize), log10(maxTargSize), nCells); % log spaced range of sizes between upper and lower test sizes in arcmin
        end
        %edited by JSk to lin space in accordance with other aim protoca
        %trialRecord(trialNo).targSize=linspace((minTargSize), (maxTargSize), nCells); % lin spaced range of sizes between upper and lower test sizes in arcmin
        
        trialRecord(trialNo).targSizePerLoc(trialRecord(trialNo).targLocs)=trialRecord(trialNo).targSize; % fill target locations with corresponding contrast
        
        trialRecord(trialNo).matchOri=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).oriErr=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).respCorrect=zeros(1,nCells); % fill all respones in each cell with zeros
        trialRecord(trialNo).stimSeen=zeros(1,nCells); % fill all seen with zeros in each cell - ignored cells counted as incorrect
        
        %imSize=maxTargSize*pixPerDeg;
        cellSizePix=pixPerDeg*max([minCellSizeDeg 1.5*maxTargSize]); % cell size - 1.5 * max target up to 6 deg
        responseRingPenWidth=pixPerDeg*0.5;
        xStimCenter=screenCenter(1)+xFactor*(cellSizePix+interCellGapPix); % x center of each cell on threen for this size target
        yStimCenter=screenCenter(2)+yFactor*(cellSizePix+interCellGapPix); % y center of each cell on threen for this size target
        
        responseRingRect=[0 0 cellSizePix-responseRingPenWidth/2 cellSizePix-responseRingPenWidth/2]; % size of the cell border
        
        for cellNum=1:nCells % draw all random orientation rings
            cSize=ceil(trialRecord(trialNo).targSizePerLoc(cellNum)*pixPerDeg);
            mylandoltC=makeLandoltC(cSize);
            mylandoltC(mylandoltC==0)=backLUT;
            mylandoltC(mylandoltC==1)=targLUT;
            trialRecord(trialNo).mylandoltC{cellNum}=mylandoltC;
            testAngle=360*rand(); % random test orientation - TODO work on this for astigmatism etc
            trialRecord( trialNo).targOri(cellNum)=testAngle; % update record for test angle
            destRect=CenterRectOnPoint([0 0 cSize cSize], xStimCenter(cellNum), yStimCenter(cellNum)); % rect for target size in this cell
            
            myTex=Screen('MakeTexture', windowPtr,mylandoltC);%, [],[],2); % write image to texture, scaled
            Screen('DrawTexture', windowPtr, myTex, [0 0 cSize cSize],destRect,testAngle); % draw texture to on screen location
            Screen('Close', myTex);
            
            destRect=CenterRectOnPoint(responseRingRect, xStimCenter(cellNum), yStimCenter(cellNum)); % add response ring
            Screen('FrameOval', windowPtr,respLUT,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
        end
        
        [~, ~, textBox] = DrawFormattedText(windowPtr,  'Click on the orientation of each C:', 100, 100, 0); % give instructions to subject
        Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw exit button
        DrawFormattedText(windowPtr, sprintf('%s', 'Next→'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
        Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
        if saveImage==1
            myImfileName=sprintf('AIMVA_LowVisionTrial%d.jpg', trialNo);
            myImage=Screen('GetImage', windowPtr);
            imwrite(myImage,myImfileName);
        end
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
            
            %             if mx < 1.5*imSize && my < 1.5*imSize % observer clicked finished area
            if mx > nextRect(1) % observer clicked finished word
                clickedExit=1;
            else
                mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                [~,respNum]=min(mouseDistFromEachBox(:)); % closest cell
                %                 mouseOriWRTChosenBox=90+atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                mouseOriWRTChosenBox=atan2d(my-yStimCenter(respNum),mx-xStimCenter(respNum)); % perceived orientation for this response
                destRect=CenterRectOnPoint(responseRingRect, xStimCenter(respNum), yStimCenter(respNum)); % rect for response feedback
                Screen('FrameOval', windowPtr,resppostIndicationLUT,destRect, responseRingPenWidth, responseRingPenWidth); % draw grey line around cells
                %                 Screen('FrameArc',windowPtr,[255 200 200],destRect,mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                Screen('FrameArc',windowPtr,respARCLUT,destRect,90+mouseOriWRTChosenBox-testGapWidthDeg/2,testGapWidthDeg,responseRingPenWidth,responseRingPenWidth); % draw orientation response arc
                Screen('Flip', windowPtr, [], 1); % show orientation and seen ring, but don't erase buffer
            end
            trialRecord( trialNo).matchOri(respNum)=mouseOriWRTChosenBox; % record reported orientation
            trialRecord( trialNo).oriErr(respNum)=diff(unwrap([trialRecord( trialNo).targOri(respNum),mouseOriWRTChosenBox]/180*pi)*180/pi); % calculate difference betwen actual and reporte orientation
            trialRecord( trialNo).stimSeen(respNum)=1; % mark this cell as response made
            %edited by JSk,2024
            %collect mouse data
            trialRecord( trialNo).MouseX=mouseX; %x coordinates for mouse
            trialRecord( trialNo).MouseY=mouseY; %y coordinates for mouse
            trialRecord( trialNo).Mousebuttons=mouseclick; %x coordinates for mouse
            
            if abs(trialRecord(trialNo).oriErr(respNum))<(360/(nAFC*2)) % if within nAFC error
                trialRecord(trialNo).respCorrect(respNum)=1; % score correct
            else
                trialRecord(trialNo).respCorrect(respNum)=0; % score incorrect
            end
        end
        Screen('Flip', windowPtr); % show stimulus and don't erase buffer
        trialRecord(trialNo).chartTime=toc(chartStart);
        trialNo=trialNo+1;
        % if saveImage==1 % saving demo image?
        %     myImfileName=sprintf('AIMAcuity%d.jpg', trialNo); % create filename
        %     myImage=Screen('GetImage', windowPtr); % grab screnshot
        %     imwrite(myImage,myImfileName); % write screenshot to image file
        % end
        totalNTrials=totalNTrials + nCells;% total # observations so far including the new charts number of cells
        
    end
    
    expDuration=toc(expStart); % stop timer - how long was experiment?
    Screen('CloseAll'); % close all windows
    
    save(dataFile,'trialRecord','expDuration','scrnWidthCm','viewDistCm','frameRate','winRect','sName','backLUT','minNumReports','timestamp_experiment_start'); % save parameters
    
    %     testLevel=[]; % set blank list of all levels for this SF
    %     respCorrect=[]; % blank list of all stimuli seen for this SF
    %     for trialSoFar=1:nTrials % work through all trials
    %         testLevel=[testLevel trialRecord(trialSoFar).targSizePerLoc]; % concatenate all test sizes
    %         respCorrect=[respCorrect  trialRecord(trialSoFar).respCorrect]; % concatenate all orientation errors
    %     end
    %     testLevel=log10((pixSizeArcmin*testLevel)/12);
    %     [acuityFitobject,~,threshAcuity]=fitPFuncnAFC_AIMLOWVISIONAcuity(testLevel, respCorrect, nAFC, 1); % use matlab's fit function to fit cumulative Gaussian to the data
    %     logMAR2Snellen=20*(10^threshAcuity); % logMAR to Snellen
    %     title([sprintf(': Participant %s Acuity: logMAR %1.2f Snellen 20/%.1f', sName, threshAcuity, logMAR2Snellen)]);
    %     xlabel('Letter Size (deg)'); % label the x axis
    %     ylabel('Proportion Correct'); % label y axis
    %     ylim([0 1]); % fix upper and lower bounds
    %     save(dataFile,'acuityFitobject','logMAR2Snellen','threshAcuity','pixSizeArcmin', 'trialRecord','expDuration','scrnWidthCm','viewDistCm','frameRate','winRect');
    
    figure();
    testLevel=[]; % set blank list of all levels for this SF
    respErr=[]; % blank list of all stimuli seen for this SF
    stimSeen=[]; % blank list of all stimuli seen for this SF
    for trialSoFar=1:nTrials % work through all trials
        testLevel=[testLevel trialRecord(trialSoFar).targSizePerLoc]; % concatenate all test sizes
        respErr=[respErr trialRecord(trialSoFar).oriErr]; % concatenate all orientation errors
        stimSeen=[stimSeen trialRecord(trialSoFar).stimSeen]; % concatenate all responses made
    end
    respErr(stimSeen==0)=ranreporterror; % assign error on ignored cells to 90 deg
    %testLevel=log10((pixSizeArcmin*testLevel)/12);%
    [angleFitobject,~,angThreshAcuity,ci]=fitPFuncOriErr_AIMLOWVISIONAcuity(testLevel',abs(respErr)',ranreporterror,1); % fit cumulative Gaussian to the ori error data
    %[angleFitobject,~,angThreshAcuity,ci]=fitPFuncOriErr(testLevel',abs(respErr)',ranreporterror,1); % fit cumulative Gaussian to the ori error data
    save(dataFile,'ExpDate','angleFitobject','angThreshAcuity', 'pixSizeArcmin', 'trialRecord','expDuration','scrnWidthCm','viewDistCm','frameRate','winRect','sName','backLUT','minNumReports','timestamp_experiment_start');
catch exception
    Screen('CloseAll');
    expDuration=toc(expStart); % stop timer - how long was experiment?
    
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'AIM_VA_LowVisionFailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sAIM_VA_LowVisionFailedRun.mat', testSName);       % file path to Mat files
    %edited Jsk, 2024
    save(dataFile,'ExpDate','exception','trialRecord','expDuration','scrnWidthCm','viewDistCm','frameRate','winRect','sName','backLUT','minNumReports','timestamp_experiment_start'); % save parameters
end


% This software is patented and owned by Northeastern University, Boston, USA; Patent number: 63/075,084 INV-21015
% Invented and developed by Peter J. Bex and Jan Skerswetat, 2021
%% Version 02/10/2022
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject ] = fitPFuncnAFC_AIMLOWVISIONAcuity( level, response, nAFC, graphData )
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
function [ fitobject, gof, thresholdEst, ci, chiSqFitReject, lettAcuityEquiv] = fitPFuncOriErr_AIMLOWVISIONAcuity( tSize, mErr,ranErr, graphData,eye,sName, criterion )
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
    set(gca, 'XScale', 'log');
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

