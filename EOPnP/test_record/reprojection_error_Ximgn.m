function [err] = reprojection_error_Ximgn(Xw,Ximgn,R,T)
% Xw ----- world coordinates of reference points by cols
% Ximgn ----- normalized image coordinates of reference points by cols
%参考点数量
n=size(Xw,2);

Xw_h=[Xw;ones(1,n)];
Ximgnrep_= [R,T] * Xw_h;

%使重投影的归一化图像坐标第三行为1
Ximgnrep = zeros(2,n);
Ximgnrep(1,:) = Ximgnrep_(1,:) ./ Ximgnrep_(3,:);
Ximgnrep(2,:) = Ximgnrep_(2,:) ./ Ximgnrep_(3,:);

%求归一化坐标的重投影误差
err_ = sqrt((Ximgn(1,:)-Ximgnrep(1,:)).^2+(Ximgn(2,:)-Ximgnrep(2,:)).^2);
err=sum(err_)/n;


