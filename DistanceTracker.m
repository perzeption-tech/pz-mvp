function distEst = DistanceTracker()
cam = webcam;
cam.Resolution = '1280x720'; %  640x480
videoFrameRGB = snapshot(cam);

standardEyeWidthAt60cm=100;
dist = 60;
% standardFaceWidthAt60cm=150;
figure('Name','Distance Tracker');
set(gcf, 'Position', get(0, 'Screensize'));

% faceDetector = vision.CascadeObjectDetector('FrontalFaceCART');
eyesDetector = vision.CascadeObjectDetector('EyePairBig');

% eye=1; view_distance='40'; Patient_Name='hvh'; Patient_ID = 1;
% f = parfeval(backgroundPool,@AIM_Acuity_SingleEye_APP, 1, eye, view_distance, Patient_Name, Patient_ID);
% disp(f.State) 
% [~, thisResult] = fetchNext(f);
for frameNo=1:4000
    videoFrameRGB = snapshot(cam);
        eyesbox = step(eyesDetector, videoFrameRGB);
        if eyesbox
            hold on
            distEst=(dist*standardEyeWidthAt60cm/eyesbox(3)) + (0.56 * (dist*standardEyeWidthAt60cm/eyesbox(3)));
            disp(distEst)
            videoFrameRGB = insertObjectAnnotation(videoFrameRGB,'rectangle',eyesbox,sprintf('Distance - %.1f cm',distEst), 'Color', 'green');
            imshow(videoFrameRGB);
            break;
        end
    
    drawnow;
    hold off
end

clear cam;
% clear all;
end

