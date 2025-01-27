function [K,W,V] = Planar_compute_K_Ximgn(Xw,Ximgn)
%compute_K  
%输入参数要求：
% Ximgn ----- 归一化图像点坐标矩阵，每列为一个点的坐标，坐标要求非齐次（2*n）
% Xw ----- 点世界坐标矩阵，每列为一个点的坐标，坐标要求非齐次（3*n）

%构造矩阵 W(2n*3) 以及矩阵 V(2n*6)

n = size(Xw,2);

w1 = repmat(eye(2),n,1);
w2 = reshape(Ximgn,2*n,1);
W = [w1 -w2];

V = zeros(2*n,6);
v0 = reshape(Xw(1:2,:)',n,2);
v3_temp1 = Ximgn(1,:).*Xw(1:2,:);
v3_temp2 = Ximgn(2,:).*Xw(1:2,:);
V(1:2:2*n-1,1:2) = -v0;
V(2:2:2*n  ,3:4) = -v0;
V(1:2:2*n-1,5:6) = v3_temp1';
V(2:2:2*n  ,5:6) = v3_temp2';


sum_ui_vi = sum(Ximgn,2);    % sum(  ,2)按行求和
sum_ui2_plus_vi2 = sum(sum(Ximgn.*Ximgn));
sum_xi_yi = sum(Xw(1:2,:),2);
sum_uixi_uiyi = sum(v3_temp1,2);
sum_vixi_viyi = sum(v3_temp2,2);
sum_ui2xi_ui2yi = sum(Ximgn(1,:).*v3_temp1,2);
sum_vi2xi_vi2yi = sum(Ximgn(2,:).*v3_temp2,2);
WTW = [ n,0,-sum_ui_vi(1); 
        0,n,-sum_ui_vi(2); 
        -sum_ui_vi(1),-sum_ui_vi(2),sum_ui2_plus_vi2];
% WTW_inv = inv(WTW);
WTV = [ -sum_xi_yi(1), -sum_xi_yi(2), 0, 0, sum_uixi_uiyi(1), sum_uixi_uiyi(2);
        0, 0, -sum_xi_yi(1), -sum_xi_yi(2), sum_vixi_viyi(1), sum_vixi_viyi(2);
        sum_uixi_uiyi(1), sum_uixi_uiyi(2), ...
        sum_vixi_viyi(1), sum_vixi_viyi(2), ...
        -( sum_ui2xi_ui2yi(1)+sum_vi2xi_vi2yi(1) ), ...
        -( sum_ui2xi_ui2yi(2)+sum_vi2xi_vi2yi(2) )];
K = W * (WTW\WTV) - V;

end

