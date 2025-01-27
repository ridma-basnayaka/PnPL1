function [R,t] = EOPnPL_ordinary(Xw,Ximgn,xs,xe,Xs,Xe)

nl = size(Xs,2);
if nl > 0
    [Rot,trans,cgrs,K,W,V] = EOPnPL_main(Xw,Ximgn,xs,xe,Xs,Xe);
    [minInd,Xsf,Xef,Lcf] = findBestRTReproj(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
    [Xsc, Xec] = moveLineCloser2DLcVector(Xsf,Xef, xs, xe, Lcf); % 输出 Xsl、Xel 是相机坐标系坐标，函数外应转换为世界坐标系坐标
    R0 = Rot(:, :, minInd);
    t0 = trans(:, minInd);    
    Xs = R0'*Xsc - repmat(R0'*t0, 1, size(Xsc, 2));
    Xe = R0'*Xec - repmat(R0'*t0, 1, size(Xsc, 2));  
end


[Rot,trans,cgrs,K,W,V] = EOPnPL_main(Xw,Ximgn,xs,xe,Xs,Xe);
[minInd,~] = Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn); %返回最小重投影误差及其索引 
R1 = Rot(:, :, minInd);
t1 = trans(:, minInd); 
s0 = cgrs(:,minInd);

s = damped_newton_EOPnPL(s0,K,W,V,Xw,Ximgn,xs, xe, Xs, Xe,15);

for i = 1 : size(s,2)
    [R(:,:,i),t(:,i),~,~] = s2Rt(s(:,i),W,V);
end

end

