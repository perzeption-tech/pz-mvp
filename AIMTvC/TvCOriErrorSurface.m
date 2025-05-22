function oriErrSurf =TvCOriErrorSurface( pedVal, cIncVal, IN, Psi, slope, minErr, guessRate, plotFigures)

% pedRange=[.01 .64];
% nCells=16;
% pedVal=logspace(log10(pedRange(1)), log10(pedRange(2)), nCells); % range of pedestal levels
% cIncVal=logspace(log10(0.002), log10(0.32), 100); % range of contrast increment levels
% [pedVal, cIncVal]=meshgrid(pedVal, cIncVal); % 2D range of contrasts and possible SFs on this system
% IN = 0.02;
% Psi = 0.5;
% slope = 0.1;
% minErr = 5;
% guessRate=90; % 90 for circles, 45 for grating
% plotFigures=1;
% 
% TvCfunc=sqrt((IN.^2+pedVal.^2).*(1+Psi))-pedVal;

% oriErrSurf=minErr + (guessRate-minErr).*(0.5-0.5.*erf((TvCfunc-cIncVal)./(sqrt(2).*slope)));

% oriErrSurf=minErr + (guessRate-minErr).*(0.5-0.5.*erf((log10(TvCfunc)-log10(cIncVal))./(sqrt(2).*slope)));
oriErrSurf=guessRate - (guessRate-minErr).*(0.5-0.5.*erf((log10(sqrt((IN.^2+pedVal.^2).*(1+Psi))-pedVal)-log10(cIncVal))./(sqrt(2).*slope)));

if plotFigures
    mesh(pedVal, cIncVal, oriErrSurf);
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlabel('log(Pedestal Contrast)');
    ylabel('log(Increment Threshold)');
    zlabel('AIM Error (deg)');
end
end
