function [K,W,V] = compute_K_Ximgn(Xw,Ximgn)
%compute_K  
%输入参数要求：
% Ximgn ----- 归一化图像点坐标矩阵，每列为一个点的坐标，坐标要求非齐次（2*n）
% Xw ----- 点世界坐标矩阵，每列为一个点的坐标，坐标要求非齐次（3*n）

%构造矩阵 W(2n*3) 以及矩阵 V(2n*9)

n = size(Xw,2);

w1 = repmat(eye(2),n,1);
w2 = reshape(Ximgn,2*n,1);
W = [w1 -w2];

V = zeros(2*n,9);
v0 = reshape(Xw',n,3);
v3_temp1 = Ximgn(1,:).*Xw;
v3_temp2 = Ximgn(2,:).*Xw;
V(1:2:2*n-1,1:3) = -v0;
V(2:2:2*n  ,4:6) = -v0;
V(1:2:2*n-1,7:9) = v3_temp1';
V(2:2:2*n  ,7:9) = v3_temp2';


sum_ui_vi = sum(Ximgn,2);    % sum(  ,2)按行求和
sum_ui2_plus_vi2 = sum(sum(Ximgn.*Ximgn));
sum_xi_yi_zi = sum(Xw,2);
sum_uixi_uiyi_uizi = sum(v3_temp1,2);
sum_vixi_viyi_vizi = sum(v3_temp2,2);
sum_ui2xi_ui2yi_ui2zi = sum(Ximgn(1,:).*v3_temp1,2);
sum_vi2xi_vi2yi_vi2zi = sum(Ximgn(2,:).*v3_temp2,2);
WTW = [ n,0,-sum_ui_vi(1); 
        0,n,-sum_ui_vi(2); 
        -sum_ui_vi(1),-sum_ui_vi(2),sum_ui2_plus_vi2];
% WTW_inv = inv(WTW);
WTV = [ -sum_xi_yi_zi(1), -sum_xi_yi_zi(2), -sum_xi_yi_zi(3), 0, 0, 0, sum_uixi_uiyi_uizi(1), sum_uixi_uiyi_uizi(2), sum_uixi_uiyi_uizi(3);
        0, 0, 0, -sum_xi_yi_zi(1), -sum_xi_yi_zi(2), -sum_xi_yi_zi(3), sum_vixi_viyi_vizi(1), sum_vixi_viyi_vizi(2), sum_vixi_viyi_vizi(3);
        sum_uixi_uiyi_uizi(1), sum_uixi_uiyi_uizi(2), sum_uixi_uiyi_uizi(3), ...
        sum_vixi_viyi_vizi(1), sum_vixi_viyi_vizi(2), sum_vixi_viyi_vizi(3), ...
        -( sum_ui2xi_ui2yi_ui2zi(1)+sum_vi2xi_vi2yi_vi2zi(1) ), ...
        -( sum_ui2xi_ui2yi_ui2zi(2)+sum_vi2xi_vi2yi_vi2zi(2) ), ...
        -( sum_ui2xi_ui2yi_ui2zi(3)+sum_vi2xi_vi2yi_vi2zi(3) )];
K = W * (WTW\WTV) - V;


end



% %构造矩阵 W(2n*3) 以及矩阵 V(2n*9)
% tic;
% n = size(Xw,2);
% W = zeros(2*n,3);
% V = zeros(2*n,9);
% for i = 0:n-1
%     W(2*i+1,1) = 1;
%     W(2*i+1,3) = -Ximgn(1,i+1);
%     W(2*i+2,2) = 1;
%     W(2*i+2,3) = -Ximgn(2,i+1);
%     
%     V(2*i+1,1:3) = -Xw(:,i+1)';
%     V(2*i+1,7:9) = Ximgn(1,i+1) .* Xw(:,i+1)';
%     V(2*i+2,4:6) = -Xw(:,i+1)';
%     V(2*i+2,7:9) = Ximgn(2,i+1) .* Xw(:,i+1)';
% end
% toc;
% %计算K
% K = W * pinv(W) * V - V;




