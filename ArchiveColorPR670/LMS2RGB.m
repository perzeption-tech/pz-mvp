% LMS2RGB.m
% This script calculates LMS cone-isolating directions.
% -------------------------------------------------------------------------

clearvars
clc

%% HP Spectra

% spectrum = csvread('LGMonitorSpectrum04-Feb-2021 11_10_51.csv');
spectrum = xlsread('Color_cali_averageing of PR670values.xlsx');
% load('SPD2019JanOct.mat');spectrum = SPD2019JanOct;
filename = 'Results_LMS2RGB_SurfacePro8_02022023.csv';

% plot
figure;
plot(spectrum(:,1),spectrum(:,2),'r','LineWidth',3);hold on
plot(spectrum(:,1),spectrum(:,3),'g','LineWidth',3);
plot(spectrum(:,1),spectrum(:,4),'b','LineWidth',3);
title('Display Spectra');legend('R','G','B');hold off

%% LMS2RGB

load SSconefund_0.mat;conefund = SSconefund_0; % SS cone fundamentals
convDirec = 2; % 1: RGB2LMS 2:LMS2RGB

inputDirec1 = [1 0 0]; % L ccVec [-1,1]
inputDirec2 = [0 1 0]; % M ccVec [-1,1]
inputDirec3 = [0 0 1]; % S ccVec [-1,1]

[outDirec1,VecProp1] = NineNumberCalc(conefund, inputDirec1, convDirec, spectrum);
[outDirec2,VecProp2] = NineNumberCalc(conefund, inputDirec2, convDirec, spectrum);
[outDirec3,VecProp3] = NineNumberCalc(conefund, inputDirec3, convDirec, spectrum);

%% Save data
L = [inputDirec1(1);inputDirec2(1);inputDirec3(1)];
M = [inputDirec1(2);inputDirec2(2);inputDirec3(2)];
S = [inputDirec1(3);inputDirec2(3);inputDirec3(3)];
VecProp = [VecProp1;VecProp2;VecProp3];

R_old = 2.*[outDirec1(1);outDirec2(1);outDirec3(1)]-1;
G_old = 2.*[outDirec1(2);outDirec2(2);outDirec3(2)]-1;
B_old = 2.*[outDirec1(3);outDirec2(3);outDirec3(3)]-1;
R_new = [outDirec1(1);outDirec2(1);outDirec3(1)];
G_new = [outDirec1(2);outDirec2(2);outDirec3(2)];
B_new = [outDirec1(3);outDirec2(3);outDirec3(3)];

T = table(L,M,S,VecProp,R_old,G_old,B_old,R_new,G_new,B_new);
writetable(T,filename);
