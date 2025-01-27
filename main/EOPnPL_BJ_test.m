clear; 
clc;
close all;

% rng('default');
addpath('util\');
addpath('..\model');
addpath('..\EOPnPL_back_project');
addpath('..\EOPnPL');

%focal length in pixels
fpix = 800;
%varying noise level
nll = 5;  nlp = 0;
%length of a line segment in 3D 
segLenBig = 3;
%scale of shift between line segments projected and present in the model
segLenShiftScale = 0.5*segLenBig;
%number of lines or points
nlines = 10;  npt = 0;
     
% camera's parameters
width= 640;
height= 480;       

[R, t, XLTw, XXw, xxn, xs,xe,Xs,Xe] = model.setup_3d_scene(npt, nlines, segLenBig, segLenShiftScale, nll, nlp, fpix);

% test EOPnPL_BJ
tic
[R1,t1] = EOPnPL_BJ(XXw,xxn,xs,xe,Xs,Xe,0,0,0);
toc
num_solutions = size(R1,3);
error_R1 = zeros(1,num_solutions);  error_t1 = zeros(1,num_solutions);
for jjj = 1:num_solutions
    tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
    error_R1(1,jjj) = tempy(1);
    error_t1(1,jjj) = tempy(2);
end


% test EOPnPL
tic
[R2,t2] = EOPnPL(XXw,xxn,xs,xe,Xs,Xe);
toc
num_solutions = size(R2,3);
error_R2 = zeros(1,num_solutions);  error_t2 = zeros(1,num_solutions);
for jjj = 1:num_solutions
    tempy = cal_pose_err([R2(:,:,jjj) t2(:,jjj)],[R t]);
    error_R2(1,jjj) = tempy(1);
    error_t2(1,jjj) = tempy(2);
end