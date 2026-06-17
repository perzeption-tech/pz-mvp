function [outDirec,veclength] = NineNumberCalc(conefund, inputDirec, convDirec, spectra)

% Inputs
%    conefund: LMS cone fundamentals; default = Stockman-Sharpe cone 
%              fundamentals (range=390~830nm, stepsize=0.1nm).
%              4 columns: wv, L, M, S
%    inputDirec: the set of directions that need to be converted, either 
%                = RGB (when convDirec=1; range[0,1]) 
%                or = LMS (when convDirec=2; range[-1,1]).
%    convDirec: =1 (RGB to LMS) =2 (LMS to RGB)
%    spectra: display spectra (range=390~780nm, stepsize=5nm)
%             4 columns: wv, R, G, B

% Outputs
%    outDirec: resulting directions, RGB or LMS
%    veclength: always cone contrast vector length


if convDirec == 1 % RGB to LMS
inputRGB = inputDirec;
elseif convDirec == 2 % LMS to RGB
inputLMS = inputDirec;
end
    
%% Energy spectra

% read in spectrum measurments
waves = spectra(:,1); % at 5 nm steps
redvals = spectra(:,2); % R
grevals = spectra(:,3); % G
bluvals = spectra(:,4); % B

% % plot data and check by eyeballing
% figure; scatter(waves,redvals,'r');hold on
% scatter(waves,grevals,'g');
% scatter(waves,bluvals,'b');hold off

%% convert energy unit to quantal unit
   % multiply by wavelength to put these into units that are proportional
   % to quanta. scaling is arbitrary (left out the physical constants).
rquantvals = waves.*redvals;
gquantvals = waves.*grevals;
bquantvals = waves.*bluvals;

% interpolation of the spectrum in quantal (discrete with 5nm step)
funrq = @(Xn) interp1(waves,rquantvals,Xn,'pchip');
fungq = @(Xn) interp1(waves,gquantvals,Xn,'pchip');
funbq = @(Xn) interp1(waves,bquantvals,Xn,'pchip');

%% Energy Fundamentals
% 
% % Stockman & Sharpe fundamentals - linear - downloaded from CVRE.ORG
% % load('energyfunds0.mat'); 
% energyfunds0 =  csvread('linss2_10e_fine.csv');
% wavele = energyfunds0(:,1); % not same range as PR650,limit scales to inner limits
% lfun = energyfunds0(:,2);
% mfun = energyfunds0(:,3);
% sfun = energyfunds0(:,4); % sfun = energyfunds0(1:226,4);
% 
% % interpolation
% funle = @(Xn) interp1(wavele,lfun,Xn,'pchip');
% funme = @(Xn) interp1(wavele,mfun,Xn,'pchip');
% funse = @(Xn) interp1(wavele,sfun,Xn,'pchip');
% 
% % wavesle2 = (390:.1:830)';
% % wavesle3 = (390:.1:615)';
% % lenergy = funle(wavesle2);
% % menergy = funme(wavesle2); 
% % senergy = funse(wavesle2); % senergy = funse(wavesle3); 

%% Quantal

% Stockman & Sharpe fundamentals - log - downloaded from CVRE.ORG
% quantalfunds0 =  csvread('ss2_10q_1.csv');
wavelq = conefund(:,1)'; % not in log, just waves
quantalfunds0 = conefund;

% since these are log, replace all the blanks with -99999.
quantalfunds0(quantalfunds0 == 0)=NaN;

lfunq = 10.^quantalfunds0(:,2); % in log, do anti-log
mfunq = 10.^quantalfunds0(:,3);
sfunq = 10.^quantalfunds0(:,4); % sfunq = 10.^quantalfunds0(1:226,4);

% interpolation
funlq = @(Xn) interp1(wavelq,lfunq,Xn,'pchip');
funmq = @(Xn) interp1(wavelq,mfunq,Xn,'pchip');
funsq = @(Xn) interp1(wavelq,sfunq,Xn,'pchip');

%% 9 numbers

% setup variables
wavestep = 5; % PR650 step
minwave = 390; % Stockman minimum wavelength
maxwave = 780 - wavestep; % stop short to avoid extrapolating

% % ENERGY %
% clear rr gg bb ll mm ss
% rr = funre;
% gg = funge;
% bb = funbe;
% 
% ll = funle;
% mm = funme;
% ss = funse;
% 
% clear rl9num gl9num bl9num rm9num gm9num bm9num rs9num gs9num bs9num
% rl9num = 0.5 * integral(@(x)rr(x).* ll(x), minwave, maxwave);
% rm9num = 0.5 * integral(@(x)rr(x).* mm(x), minwave, maxwave);
% rs9num = 0.5 * integral(@(x)rr(x).* ss(x), minwave, maxwave);
% gl9num = 0.5 * integral(@(x)gg(x).* ll(x), minwave, maxwave);
% gm9num = 0.5 * integral(@(x)gg(x).* mm(x), minwave, maxwave);
% gs9num = 0.5 * integral(@(x)gg(x).* ss(x), minwave, maxwave);
% bl9num = 0.5 * integral(@(x)bb(x).* ll(x), minwave, maxwave);
% bm9num = 0.5 * integral(@(x)bb(x).* mm(x), minwave, maxwave);
% bs9num = 0.5 * integral(@(x)bb(x).* ss(x), minwave, maxwave);
% 
% nineNumbsEnergy = [rl9num,gl9num,bl9num;rm9num,gm9num,bm9num;rs9num,gs9num,bs9num];

% Quantal %
clear rr gg bb ll mm ss
rr = funrq;
gg = fungq;
bb = funbq;

ll = funlq;
mm = funmq;
ss = funsq;

clear rl9num gl9num bl9num rm9num gm9num bm9num rs9num gs9num bs9num
rl9num = 0.5 * integral(@(x)rr(x).* ll(x), minwave, maxwave);
rm9num = 0.5 * integral(@(x)rr(x).* mm(x), minwave, maxwave);
rs9num = 0.5 * integral(@(x)rr(x).* ss(x), minwave, maxwave);
gl9num = 0.5 * integral(@(x)gg(x).* ll(x), minwave, maxwave);
gm9num = 0.5 * integral(@(x)gg(x).* mm(x), minwave, maxwave);
gs9num = 0.5 * integral(@(x)gg(x).* ss(x), minwave, maxwave);
bl9num = 0.5 * integral(@(x)bb(x).* ll(x), minwave, maxwave);
bm9num = 0.5 * integral(@(x)bb(x).* mm(x), minwave, maxwave);
bs9num = 0.5 * integral(@(x)bb(x).* ss(x), minwave, maxwave);

nineNumbsQuantal = [rl9num,gl9num,bl9num;rm9num,gm9num,bm9num;rs9num,gs9num,bs9num];
% save('nineNumbsQuantal','nineNumbsQuantal')
% Check %
% check = nineNumbsQuantal./nineNumbsEnergy;

%% Using the 9number matrix
   % We want to go from an rgbVector, consisting of relative modulations of 
   % guns +- {1,1,1}, to cone contrasts. An rgbVector of {0,0,0} represents
   % the mean state - grey, midluminance. So does the corresponding ccVector
   % of {0,0,0}. Want to define the rgbVectors, not on the unit sphere,
   % but the unit cube so we aren't wasting energy we could use. We also want
   % to go from ccVector to rgbVector. Normally we'll do that in relative terms,
   % to get the relative rgbVector for a particular cone contrast direction,
   % and use it by scaling by contrast (plus or minus).

% find mean vector (grey), unit meaningful
meanLMS = (nineNumbsQuantal*[1 1 1]')';
    % the corresponding meanRGB is [0 0 0].

if convDirec == 1    %%% go from rgbVector to ccVector %%%
% INPUT new RGB [0,1]
% convert to old RGB[-1,1]
rContrOLD = (inputRGB(1)*2)-1; 
gContrOLD = (inputRGB(2)*2)-1;
bContrOLD = (inputRGB(3)*2)-1;
rgbVec = [rContrOLD gContrOLD bContrOLD]; % real RGB contrast

ccVec = (nineNumbsQuantal*rgbVec')'./meanLMS; % [GET THIS]
 % function to calculate triplet for a given RGB vec.
 % results explanation: + means increase, - means decrese, number indicates
 % magnitude. ex. [-0.12;-0.19;-0.90], meaning that this RGB (relative to mean)
 % produces a 12% decrement in L-cones, 19% decrement in M-cones, 90%
 % decremnet in S-cones.
 
ccVecNorm = ccVec/max(abs(ccVec)); % [GET THIS] %LMS[-1,1] 
ccVecLength = ((ccVec.^2)*[1 1 1]')^0.5; % [GET THIS]

outDirec = ccVecNorm; % Function output 1
veclength = ccVecLength; % Function output 2

% filename_vec = strcat('VecLength_',string(outDirec(1)),'.mat');
% save(filename_vec,'veclength');

elseif convDirec == 2  %%% go from ccVector to rgbVector %%%
% INPUT
ccVec2 = [inputLMS(1) inputLMS(2) inputLMS(3)]; % LMS direction [-1,1] 
rgbVec2 = (nineNumbsQuantal^(-1)*(meanLMS.*ccVec2)')'; 
 % (rgb vector length defaults to 1, but could be different)
rgbVecNorm = rgbVec2/max(abs(rgbVec2)); % old RGB [-1,1];normalize to unit cube
rgbVecNormnew = (rgbVecNorm.*.5)+.5; % new RGB [0,1]
% rgbVecLength = ((ccVec2.^2)*[1 1 1]')^0.5;%rgbVecLength=1

outDirec = rgbVecNormnew; % Function output 1

% get CC vec length
ccVec = (nineNumbsQuantal*rgbVecNorm')'./meanLMS; % [GET THIS]
 % results explanation: + means increase, - means decrese, number indicates
 % magnitude. if [-0.12;-0.19;-0.90], meaning that this RGB (relative to mean)
 % produces a 12% decrement in L-cones, 19% decrement in M-cones, 90%
 % decremnet in S-cones.
 
ccVecNorm = ccVec/max(abs(ccVec)); % LMS[-1,1] 
ccVecLength = ((ccVec.^2)*[1 1 1]')^0.5; % [GET THIS]
veclength = ccVecLength;  % Function output 2

end

end