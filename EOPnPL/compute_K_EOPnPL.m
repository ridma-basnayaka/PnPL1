function [K,W,V] = compute_K_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe)
% 输入参数要求
% Xwp: 参考点 世界坐标系坐标，非齐次，每列一点
% Ximgnp: 参考点 归一化图像像素坐标，非齐次，每列一点
% Xs: 线起始点 世界坐标系坐标，非齐次，每列一点
% Xe: 线结束点 世界坐标系坐标，非齐次，每列一点
% xs: 线起始点 归一化图像像素坐标，非齐次，每列一点
% xe: 线结束点 归一化图像像素坐标，非齐次，每列一点


% 只有线特征
if (isempty(Xw) && isempty(Ximgn))
    
%     nl = size(Xs,2);
%     Wl = zeros(2*nl,3);
%     Vl = zeros(2*nl,9);
%     I = pnl_preprocess(xs,xe);  
%     Pi = zeros(3,1);
%     Qi = zeros(3,1);
%     for i = 1 : nl
%         Pi = Xs(:,i);
%         Qi = Xe(:,i);
%         Wl( 2*i-1 : 2*i , : ) = repmat(I(:,i)',2,1);
%         Vl( 2*i-1 : 2*i , 1 : 3) = -I(1,i) * [Pi';Qi'];
%         Vl( 2*i-1 : 2*i , 4 : 6) = -I(2,i) * [Pi';Qi'];
%         Vl( 2*i-1 : 2*i , 7 : 9) = -I(3,i) * [Pi';Qi'];    
%     end
% 
%     W = Wl;
%     V = Vl;
%     K = W * pinv(W) * V - V;
   

    n = size(Xs,2);
    W = zeros(2*n,3);
    V = zeros(2*n,9);
    I = pnl_preprocess(xs,xe)';   
    W(1:2:2*n-1,:) = I;
    W(2:2:2*n  ,:) = I;
    V(1:2:2*n-1,:) = -[I(:,1).*(Xs') I(:,2).*(Xs') I(:,3).*(Xs')];
    V(2:2:2*n  ,:) = -[I(:,1).*(Xe') I(:,2).*(Xe') I(:,3).*(Xe')];
        
    WTW_11_12_13 = sum(W(:,1).*W);
    WTW_22_23 = sum(W(:,2).*W(:,2:3));
    WTW_33 = sum(W(:,3).*W(:,3));    
    WTW = [ WTW_11_12_13 ;  WTW_11_12_13(2) WTW_22_23; WTW_11_12_13(3) WTW_22_23(2) WTW_33];
%     WTW_inv = inv(WTW);
    
    WTV_11_to_19 = sum(W(:,1).*V);
    WTV_24_to_29 = sum(W(:,2).*V(:,4:9));
    WTV_37_to_39 = sum(W(:,3).*V(:,7:9));
    WTV = [ WTV_11_to_19;...
            WTV_11_to_19(4:6), WTV_24_to_29;...
            WTV_11_to_19(7:9), WTV_24_to_29(4:6), WTV_37_to_39];
    
    K = W * (WTW\WTV) - V;


% 只有点特征
elseif (isempty(Xs) && isempty(Xe))
%     np = size(Xw,2);
%     Wp = zeros(2*np,3);
%     Vp = zeros(2*np,9);
%     for i = 1 : np
%         Wp(2*i-1,1) = 1;
%         Wp(2*i-1,3) = -Ximgn(1,i);
%         Wp(2*i,2) = 1;
%         Wp(2*i,3) = -Ximgn(2,i);
% 
%         Vp(2*i-1,1:3) = -Xw(:,i)';
%         Vp(2*i-1,7:9) = Ximgn(1,i) .* Xw(:,i)';
%         Vp(2*i,4:6) = -Xw(:,i)';
%         Vp(2*i,7:9) = Ximgn(2,i) .* Xw(:,i)';
%     end
%     W = Wp;
%     V = Vp;
%     K = W * pinv(W) * V - V;
  

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
%     WTW_inv = inv(WTW);
    WTV = [ -sum_xi_yi_zi(1), -sum_xi_yi_zi(2), -sum_xi_yi_zi(3), 0, 0, 0, sum_uixi_uiyi_uizi(1), sum_uixi_uiyi_uizi(2), sum_uixi_uiyi_uizi(3);
            0, 0, 0, -sum_xi_yi_zi(1), -sum_xi_yi_zi(2), -sum_xi_yi_zi(3), sum_vixi_viyi_vizi(1), sum_vixi_viyi_vizi(2), sum_vixi_viyi_vizi(3);
            sum_uixi_uiyi_uizi(1), sum_uixi_uiyi_uizi(2), sum_uixi_uiyi_uizi(3), ...
            sum_vixi_viyi_vizi(1), sum_vixi_viyi_vizi(2), sum_vixi_viyi_vizi(3), ...
            -( sum_ui2xi_ui2yi_ui2zi(1)+sum_vi2xi_vi2yi_vi2zi(1) ), ...
            -( sum_ui2xi_ui2yi_ui2zi(2)+sum_vi2xi_vi2yi_vi2zi(2) ), ...
            -( sum_ui2xi_ui2yi_ui2zi(3)+sum_vi2xi_vi2yi_vi2zi(3) )];
    K = W * (WTW\WTV) - V;

    
    
% 点线特征均有    
% 构造点特征矩阵 Wp Vp
%构造矩阵 Wp(2np*3) 以及矩阵 Vp(2np*9)    
else
%     np = size(Xw,2);
%     Wp = zeros(2*np,3);
%     Vp = zeros(2*np,9);
%     for i = 1 : np
%         Wp(2*i-1,1) = 1;
%         Wp(2*i-1,3) = -Ximgn(1,i);
%         Wp(2*i,2) = 1;
%         Wp(2*i,3) = -Ximgn(2,i);
% 
%         Vp(2*i-1,1:3) = -Xw(:,i)';
%         Vp(2*i-1,7:9) = Ximgn(1,i) .* Xw(:,i)';
%         Vp(2*i,4:6) = -Xw(:,i)';
%         Vp(2*i,7:9) = Ximgn(2,i) .* Xw(:,i)';
%     end
% 
% % 构造线特征矩阵 Wl Vl
% %构造矩阵 Wl(2nl*3) 以及矩阵 Vl(2nl*9)
%     nl = size(Xs,2);
%     Wl = zeros(2*nl,3);
%     Vl = zeros(2*nl,9);
%     I = pnl_preprocess(xs,xe);  
%     Pi = zeros(3,1);
%     Qi = zeros(3,1);
%     for i = 1 : nl
%         Pi = Xs(:,i);
%         Qi = Xe(:,i);
%         Wl( 2*i-1 : 2*i , : ) = repmat(I(:,i)',2,1);
%         Vl( 2*i-1 : 2*i , 1 : 3) = -I(1,i) * [Pi';Qi'];
%         Vl( 2*i-1 : 2*i , 4 : 6) = -I(2,i) * [Pi';Qi'];
%         Vl( 2*i-1 : 2*i , 7 : 9) = -I(3,i) * [Pi';Qi'];    
%     end
% 
%     W = [Wp;Wl];
%     V = [Vp;Vl];
%     K = W * pinv(W) * V - V;

    
%点特征      
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
%     WTW_inv = inv(WTW);
    WTV = [ -sum_xi_yi_zi(1), -sum_xi_yi_zi(2), -sum_xi_yi_zi(3), 0, 0, 0, sum_uixi_uiyi_uizi(1), sum_uixi_uiyi_uizi(2), sum_uixi_uiyi_uizi(3);
            0, 0, 0, -sum_xi_yi_zi(1), -sum_xi_yi_zi(2), -sum_xi_yi_zi(3), sum_vixi_viyi_vizi(1), sum_vixi_viyi_vizi(2), sum_vixi_viyi_vizi(3);
            sum_uixi_uiyi_uizi(1), sum_uixi_uiyi_uizi(2), sum_uixi_uiyi_uizi(3), ...
            sum_vixi_viyi_vizi(1), sum_vixi_viyi_vizi(2), sum_vixi_viyi_vizi(3), ...
            -( sum_ui2xi_ui2yi_ui2zi(1)+sum_vi2xi_vi2yi_vi2zi(1) ), ...
            -( sum_ui2xi_ui2yi_ui2zi(2)+sum_vi2xi_vi2yi_vi2zi(2) ), ...
            -( sum_ui2xi_ui2yi_ui2zi(3)+sum_vi2xi_vi2yi_vi2zi(3) )];
     Wp = W;
     Vp = V;
     WTWp = WTW;
     WTVp = WTV;

    %线特征
    n = size(Xs,2);
    W = zeros(2*n,3);
    V = zeros(2*n,9);
    I = pnl_preprocess(xs,xe)';   
    W(1:2:2*n-1,:) = I;
    W(2:2:2*n  ,:) = I;
    V(1:2:2*n-1,:) = -[I(:,1).*(Xs') I(:,2).*(Xs') I(:,3).*(Xs')];
    V(2:2:2*n  ,:) = -[I(:,1).*(Xe') I(:,2).*(Xe') I(:,3).*(Xe')];

    WTW_11_12_13 = sum(W(:,1).*W);
    WTW_22_23 = sum(W(:,2).*W(:,2:3));
    WTW_33 = sum(W(:,3).*W(:,3));    
    WTW = [ WTW_11_12_13 ;  WTW_11_12_13(2) WTW_22_23; WTW_11_12_13(3) WTW_22_23(2) WTW_33];
    %     WTW_inv = inv(WTW);

    WTV_11_to_19 = sum(W(:,1).*V);
    WTV_24_to_29 = sum(W(:,2).*V(:,4:9));
    WTV_37_to_39 = sum(W(:,3).*V(:,7:9));
    WTV = [ WTV_11_to_19;...
            WTV_11_to_19(4:6), WTV_24_to_29;...
            WTV_11_to_19(7:9), WTV_24_to_29(4:6), WTV_37_to_39];

     Wl = W;
     Vl = V;
     WTWl = WTW;
     WTVl = WTV;

     W = [Wp;Wl];
     V = [Vp;Vl];
     WTW = WTWp + WTWl;
     WTV = WTVp + WTVl;
     K = W * (WTW\WTV) - V;

end

end

