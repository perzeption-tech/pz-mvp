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

