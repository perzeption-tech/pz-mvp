    function [fitobject] =fitTvCOriErrorSurface(PedVal, CSVal, OriErrors)

% myFitModel = fittype( 'TvCOriErrorSurface( x, y, IN, Psi, slope, minErr, 45, 0)',...
myFitModel = fittype( '45 - (45-minErr).*(0.5-0.5.*erf((log10(sqrt((IN.^2+x.^2).*(1+Psi))-x)-log10(y))./(sqrt(2).*slope)))',...
     'independent', {'x','y'}, ...
     'coefficients',{'IN','Psi','slope','minErr'},...
     'dependent','z');

fitobject = fit( [PedVal(:), CSVal(:)], OriErrors(:), myFitModel, 'StartPoint', [0.02 0.5 0.1 5], 'Lower',[0.01 0.01 0.01 1], 'Upper',[0.1 2 0.2 20]);% , 'Weights', fitWeights, 'Lower',[min(xVals) min(yObs) 0.1],'Upper',[max(xVals) max(yObs) 0.75*(max(xVals)-min(xVals))]);

end




