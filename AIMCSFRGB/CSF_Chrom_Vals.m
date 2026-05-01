if LMSType==1  % luminance defined
                        rgbGabor(:,:,1)=targContrast*myTarg; %L-M
                        rgbGabor(:,:,2)=targContrast*myTarg; %L-M
                        rgbGabor(:,:,3)=targContrast*myTarg; %L-M
                     elseif LMSType==2 % red/green modulation
                        lmsContrast=sqrt((targContrast.^2)/2);
                        %   rgbGabor(:,:,1)=lmsContrast*(rgbRatios(1,1)-rgbRatios(2,1))*myTarg; %L-M
                        % rgbGabor(:,:,2)=lmsContrast*(rgbRatios(1,2)-rgbRatios(2,2))*myTarg; %L-M
                        % rgbGabor(:,:,3)=lmsContrast*(rgbRatios(1,3)-rgbRatios(2,3))*myTarg; %L-M
%updated by JSk, April 2026
                        rgbGabor(:,:,1)=lmsContrast*(rgbRatios(3,1)-rgbRatios(4,1))*myTarg; %L-M
                        rgbGabor(:,:,2)=lmsContrast*(rgbRatios(3,2)-rgbRatios(4,2))*myTarg; %L-M
                        rgbGabor(:,:,3)=lmsContrast*(rgbRatios(3,3)-rgbRatios(4,3))*myTarg; %L-M
                        % cone contrast = sqrt(L^2 + M^2 + S^2); Kim et al (2018). A Normative Data Set for the Clinical Assessment of Achromatic and Chromatic Contrast Sensitivity Using a qCSF Approach. Investigative Opthalmology & Visual Science, 58(9), 3628.
                    elseif LMSType==3 % blue yellow modulation
                        % Cc^2 = 3*Clms^2; (where L=M=S=1)
                        % (Cc^2)/2 = Clms^2;
                        % Clms=sqrt((Cc^2)/2)
                        lmsContrast=sqrt((targContrast.^2)/3);
                        rgbGabor(:,:,1)=lmsContrast*((rgbRatios(3,1)+rgbRatios(4,1))/2-rgbRatios(2,1))*myTarg; %(L+M)/2-S
                        rgbGabor(:,:,2)=lmsContrast*((rgbRatios(3,2)+rgbRatios(4,2))/2-rgbRatios(2,2))*myTarg; %(L+M)/2-S
                        rgbGabor(:,:,3)=lmsContrast*((rgbRatios(3,3)+rgbRatios(4,3))/2-rgbRatios(2,3))*myTarg; %(L+M)/2-S
                        % Updated by JSk, April 2026
                    elseif LMSType==4 % L+ L- modulation
                        lmsContrast=sqrt((targContrast.^2)/3);
                        rgbGabor(:,:,1)=lmsContrast*(rgbRatios(3,1))*myTarg; %L
                        rgbGabor(:,:,2)=lmsContrast*(rgbRatios(3,2))*myTarg; %L
                        rgbGabor(:,:,3)=lmsContrast*(rgbRatios(3,3))*myTarg; %L
                    elseif LMSType==5 % M+ M- modulation
                                                lmsContrast=sqrt((targContrast.^2)/3);
                        rgbGabor(:,:,1)=lmsContrast*(rgbRatios(4,1))*myTarg; %M
                        rgbGabor(:,:,2)=lmsContrast*(rgbRatios(4,2))*myTarg; %M
                        rgbGabor(:,:,3)=lmsContrast*(rgbRatios(4,3))*myTarg; %M
                end