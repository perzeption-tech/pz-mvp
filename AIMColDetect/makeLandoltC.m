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

