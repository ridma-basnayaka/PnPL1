% clear
addpath others;

%参考点数量
npt = 100;
nl = 2;

% camera's parameters
width= 640;
height= 480;
f= 800;

% generate 3d coordinates in camera space
Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
t= mean(Xc,2);
R= rodrigues(randn(3,1));
XXw= inv(R)*(Xc-repmat(t,1,npt));

% projection
xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f; %相机坐标系中的投影方程,不考虑主点偏置
% xx= [f.*(Xc(1,:)./Xc(3,:)) + width/2 ; f.*(Xc(2,:)./Xc(3,:)) + height/2]; %相机坐标系中的投影方程,考虑主点偏置
% A = [f 0 width/2; 0 f height/2; 0 0 1];   %考虑主点偏置

xxn= xx+randn(2,npt)*nl;    %添加噪声
t1 = tic;
[Re,te] = EOPnP(XXw,xxn/f);     %不考虑主点偏置
% [Re,te] = EOPnP(XXw,inv(A)*[xxn;ones(1,npt)]);  %考虑主点偏置
t2 = toc(t1);
error = inf;
for jjj = 1:size(Re,3)
    tempy = cal_pose_err([Re(:,:,jjj),te(:,jjj)],[R t]);
    if sum(tempy) < error
        y = tempy;
        error = sum(tempy);
    end
end
y
