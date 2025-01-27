function [R,t] = EOPnPL_BJ_ordinary_5(Xw,Ximgn,xs,xe,Xs,Xe,iter1,iter2,PnP_iter)


nlines = size(Xs,2);

if nlines > 0
    [Rot, trans,cgrs,K,W,V] = EOPnPL_BJ_main_length_weight(Xw,Ximgn,xs,xe,Xs,Xe);
    [minInd,Xsf,Xef,Lcf] = EOPnPL_BJ_findBestRTReproj(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
    R_temp = Rot(:, :, minInd);  t_temp = trans(:, minInd);
    if iter1 ~= 0
        s_temp = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe,iter1);
        [R_temp,t_temp] = EOPnPL_BJ_s2Rt(s_temp,W,V);
    end
    [Xsc, Xec] = EOPnPL_BJ_moveLineCloser2DLcVector(Xsf,Xef, xs, xe, Lcf); % 输出 Xsl、Xel 是相机坐标系坐标，函数外应转换为世界坐标系坐标
    Xs = R_temp'*Xsc - repmat(R_temp'*t_temp, 1, size(Xsc, 2));
    Xe = R_temp'*Xec - repmat(R_temp'*t_temp, 1, size(Xsc, 2));  
end

[Rot,trans,cgrs,K,W,V] = EOPnPL_BJ_main_length_weight(Xw,Ximgn,xs,xe,Xs,Xe);
[minInd,~] = EOPnPL_BJ_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn); %返回最小重投影误差及其索引 
s = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe, 15);[R,t] = s2Rt(s,W,V);


end






