clear
addpath others;

%参考点数量
npt = 100;
nl = 2;

% camera's parameters
width= 640;
height= 480;
f= 800;

% generate 3d coordinates in camera space
XXw= [xrand(2,npt,[-2 2]); zeros(1,npt)];
R= rodrigues(randn(3,1));
% disp(R);
t= [rand-0.5;rand-0.5;rand*8+4];
% disp(t);
Xc= R*XXw+repmat(t,1,npt);

% projection
xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
xxn= xx+randn(2,npt)*nl;

tic;
[Re,te] = EOPnP(XXw,xxn/f);
toc;

error = inf;
for jjj = 1:size(Re,3)
    tempy = cal_pose_err([Re(:,:,jjj),te(:,jjj)],[R t]);
    if sum(tempy) < error
        y = tempy;
        error = sum(tempy);
    end
end
y
