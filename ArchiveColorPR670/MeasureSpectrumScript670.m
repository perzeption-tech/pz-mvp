% Script to measure spectral power distribution using the PR-655
% and save the values in a text file.
%
% NOTE: 
%       10/25/2016rte  Using BRAINARD'S convention (for ease of use)
%                       he multiplies the raw energy values by 4
%                       (so they are per wavelength BAND, not per
%                       wavelength).  Then he interpolates them to 
%                       5 nm.  CHANGE LATER

% wavelengths=linspace(380, 380+5*81, 81);
% [spd, qual] = PR650measspd(S,syncMode);  ONCE INITIALIZED, CAN JUST CALL
% THIS DIRECTLY if don't care about other stuff.

Screen('Preference', 'SkipSyncTests', 1);
sca;
close all;
clearvars;
commandwindow;

infoinput = inputdlg({'Measurement Number ='},':)',[1 50],{'0'});
nMeasure = cell2mat(infoinput);

% filename='SpectrumMeasurement1_HP';       % CHANGE THIS STRING
% filename=strcat(filepath, filename);  % specifying file path is not
% working; save to Current Folder than MOVE the file to shared folder

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

%  here setup  some globals for the experiment
global EXPWINDOW


midGrey=[.5 .5 .5];  % needed here to open background
% oldPriority=Priority(9);
% try


screenNumber=max(Screen('Screens'));
% Setup imaging pipeline:
PsychImaging('PrepareConfiguration');

% % see if using BITS#  here and prepare it if so:
% usingBitsSharpFlag = false;

% if usingBitsSharpFlag == true   
% OpenQ=SetupBitsSharp(1,'/dev/cu.usbmodem72124141');  % this string may need to change, or possibly isn't needed at all
%                         % note mode of 1 causes aspect ratio to be correct,
%                         % BUT even columns of stimulus do not appear!
% %                         OpenQ=SetupBitsSharp(1);
% end  % calling SetupBitsSharp

% Name Monitor here to help avoid confusion!
MonitorName='SurfacePro8_02022023';

usingPR670Flag = true;   % change to true if really measuring and it is connected

% pr650 actually returns 101 steps from 380-780 with bandwidth of 4;
 nWaves=81;   % using brainard's routines unmodified
 wavelengths=linspace(380, 380+5*(nWaves-1),nWaves);

%  %  turn on PR650 no more than 5 secs before initialization, otherwise
%  it will revert to manual mode
if usingPR670Flag == true
     CMCheckInit(5) ;
     % this initializes the serial port
     % meterType 4 is PR655
     % meterType 5 is PR670
end


% %%%%%%%%%%%%%%%%%%%%%%

% 
adaptSecs=1; % in seconds; this can be used for initial wait delay before starting PR650 measurements

blankPix=300;         %  rectangular area where the test will appear defined by blankPix

% %%%%%%%%%%%%%%%%%%%%%%

% 
% AssertOpenGL;
% 
% %Set higher DebugLevel, so that you don't get all kinds of messages flashed
%     %at you each time you start the experiment:
%     olddebuglevel=Screen('Preference', 'VisualDebuglevel', 3);
%     %rte: this sets the new screen preference to VisualDebuglevel=3, and
%     %saves the original preference setting in olddebuglevel so it can be
%     %restored at the end of this script.
% 
% % %%%%%%%%%%%%%%%%%%%%%%

AssertOpenGL;

%Set higher DebugLevel, so that you don't get all kinds of messages flashed
    %at you each time you start the experiment:
    olddebuglevel=Screen('Preference', 'VisualDebuglevel', 3);
    %rte: this sets the new screen preference to VisualDebuglevel=3, and
    %saves the original preference setting in olddebuglevel so it can be
    %restored at the end of this script.

% prepare output matrix:
nCols=1+2*4; % number of spd's will measure is 3 NO 4 incl white, quality column, plus one wavelegnth column
outputData=zeros(nWaves, nCols);  % will hold wavelength, red, redqual, green, green qual, blue, blue qual,in that order
outputData(:,1)=wavelengths';

try
% Enable unified mode of KbName, so KbName accepts identical key names on
%     % all operating systems (not absolutely necessary, but good practice):
%     KbName('UnifyKeyNames');  % shouldn't be necessary since did PsychDefaultSetup(2)
% 
%    
%     KbCheck;
    
screenSize = get(0,'screensize');
%find screens
    screenNumber=max(Screen('Screens'));
    [EXPWINDOW, expWindowRect] = PsychImaging('OpenWindow', screenNumber, midGrey);  % mid gray but Gamma correction!
    interframeInterval=Screen('GetFlipInterval', EXPWINDOW);  % will use for timing purposes
 
    Screen('BlendFunction', EXPWINDOW, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  % do this iff doing it in exper script


    HideCursor 




 %disable output of keypresses to Matlab. !!!use with care!!!!!!
    %if the program gets stuck you might end up with a dead keyboard
    %if this happens, press CTRL-C to reenable keyboard handling -- it is
    %the only key still recognized.
    ListenChar(2);
    
    
    %Speak('Press Y');
%       advanceKey = KbName('Y');  % keypress before the measurements begin
%      while 1
%          [a,b,keyCode] = KbCheck;
%          if keyCode(advanceKey)
%              break;
%  
%          end
%   
%        end


WaitSecs(adaptSecs) ; % set this long engough to, for example, leave room, turn off lights, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is where we will load in our gamma table\          % IMPORTANT NOTE:
%                        may need to try different sized normalized gamma
%                        tables!  using the biggest possible here
%                        now...rte10/14/16
[originalGammaClut,dacbits,reallutsize]=Screen('ReadNormalizedGammaTable',screenNumber)
Screen('LoadNormalizedGammaTable', EXPWINDOW, linspace(0, 1, 2^dacbits-1)' * [1, 1, 1], 0);        % try an expanded LUT?????
% dacbits is the number of bits of intensity resolution of the screen
% (nominally)

% here is where we would read in the gamma table IF we're doing a gamma check!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%% for debugging set these here to constant so don't need to hook up
% the PR650.  
if usingPR670Flag == false
    spd=-rand(nWaves,1);        % make'em negative so its is clear this isn't real data!
    qual=-randi(100,1);
end
% %%%%%%%%%%

% testDirection=[ 0 0 0] % just for debugging

patchMat = zeros(blankPix, blankPix, 3); % 

% now will have one loops:
% outer loop varies primary/gun (red, green, blue)


for gun = 1:4  % modulate guns
    testDirection=[0 0 0]; % reset it  FOR REAL MEASUREMENT MAKE THIS [0 0 0] ; could also do [1 1 1] and variations to test for additivity
    if gun<4
        testDirection(gun)=1;  % set this to the desired gun
    else
        testDirection=[1 1 1];
    end
    
    
        
        stimRGB=testDirection;

        Beeper; 

     
     
     %   use  FillRect iff _experiment_ uses FillRect or similar.  Use Textures if
     %   _experiment_ uses texture patches
       Screen('FillRect', EXPWINDOW, stimRGB, CenterRect([0 0 blankPix blankPix], expWindowRect));

     
%         patchMat=StraightEdgeMatrixForBoynton(stimRGB, stimRGB, 0,blankPix, blankPix, 0);  % use as much of exper routines as possible
%         patchTexture=Screen('MakeTexture',EXPWINDOW,patchMat, [],[],1);  % the final '1' is for high-resolution intensities
%        Screen('DrawTexture', EXPWINDOW, patchTexture)
        presentTime=Screen('Flip', EXPWINDOW);  %  start time  -- not needed here
        WaitSecs(1);
%         %%%%%%%%%%%%%%% PR650 READING DONE HERE
        % DON'T USE THIS SECTION IFF NOT USING PR650 (ie, for debugging )
        if usingPR670Flag==true        %
        [spd, qual] = PR670measspd;  % values are in COLUMNS (?????)
        qual=repmat(qual,nWaves,1);  % make a set of (identical) quality code values so have them...
        % spd=spd';
        else
            spd=gun*1.5*spd;  % shift vertically so can see them in plot, for debugging
            qual=gun*1.5*repmat(qual,nWaves,1);
        end
        
        
        
        %%%%%%%%%%%%%%% PR650 READING DONE HERE

       outputData(:,2*gun)=spd;
       outputData(:,1+2*gun)=qual;
       qual=-randi(100,1);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %  
       %  no need to wait around if using PR650

%         if usingPR655Flag==true
%                 Speak('Press y')
%                  while 1                 % wait here until 'Y' pressed
%                  [a,b,keyCode] = KbCheck;
%                      if keyCode(advancekey)
%                          break;
%                      end
%                  end
%         end  % end if for pr650 flag
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         Screen('Close', patchTexture);          % I think this should be done...?
    
    
end % gun loop
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FlushEvents('keyDown');	% discard all the chars from the Event Manager queue.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output the readings to a text file

whenStr=datestr(now,2);
OutFileName = strcat(MonitorName,'Spectrum',nMeasure,'.txt');
            % create the new file name as a concantenation
            % of the string 'Spectrum' and the current date and time
fileId = fopen(OutFileName, 'w');  % Open file to write
tabsOut = repmat('%f\t',1,(nCols-1)); % will put a tab after each floating
                                      % point except last
% fprintf(fid,[tabsOut,'%s\n'],colHeaders)        % first put in the header
% text  FIX THIS LATER

for row=1:size(outputData,1)
    fprintf(fileId,[tabsOut,'%f\n'],outputData(row,:));  % use fprintf to write data file
end  % the \n is a "newline" or return
fclose(fileId);

% %%%%%%%%%%%%%%%%%%%%%%%

%  plot the measurements
xx=outputData(:,1);
red=outputData(1:nWaves, 2);  % first set of meas is red
green=outputData(:, 4);  % 2 set of meas is green  % skipping qual code
blue=outputData(:,6);  % 3 set of meas is blue
figure
plot(xx,red,'r', xx,green,'g', xx, blue,'b')

    ShowCursor;
   %   sca; %or Screen('CloseAll');
    ListenChar(0);
    %return to olddebuglevel
    Screen('Preference', 'VisualDebuglevel', olddebuglevel);
    Screen('LoadNormalizedGammaTable', screenNumber, originalGammaClut);  % put it back the way it was
    
    Screen('FillRect', EXPWINDOW, [0 0 0], expWindowRect);  % black the screen and wait here until Y is pressed
    Screen('Flip', EXPWINDOW);
    
     sca;

catch exception
    % This section is executed only in case an error happens in the
    % experiment code implemented between try and catch...
    ShowCursor;
    Screen('CloseAll'); %or sca
    ListenChar(0);
    Screen('Preference', 'VisualDebuglevel', olddebuglevel);
    Screen('LoadNormalizedGammaTable', screenNumber, originalGammaClut);  % put it back the way it was
    %output the error message
    psychrethrow(psychlasterror);
end
