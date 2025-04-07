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

if plotFigures
    mesh(logSF, logCS, oriErrSurf);
    xlabel('log Spatial Frequency');
    ylabel('log Contrast');
    zlabel('AIM Error (deg)');
end
end
