function [R,t] = EOPnPL_Pm_ordinary(Xw,Ximgn,xs,xe,Xs,Xe)

        np = size(Xw,2);
        nl = size(Xs,2);
        if nl > 0
            [Rot,trans] = EOPnPL_Pm_main(Xw,Ximgn,xs,xe,Xs,Xe);
            [minInd,Xsf,Xef,Lcf] = EOPnPL_Pm_findBestRTReproj(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
            [Xsc, Xec] = EOPnPL_Pm_moveLineCloser2DLcVector(Xsf,Xef, xs, xe, Lcf); % 输出 Xsl、Xel 是相机坐标系坐标，函数外应转换为世界坐标系坐标
            R0 = Rot(:, :, minInd);
            t0 = trans(:, minInd);    
            Xs = R0'*Xsc - repmat(R0'*t0, 1, size(Xsc, 2));
            Xe = R0'*Xec - repmat(R0'*t0, 1, size(Xsc, 2));  
        end

        
        [Rot,trans,cgrs,K,W,V] = EOPnPL_Pm_main(Xw,Ximgn,xs,xe,Xs,Xe);
        if  ( np + nl ) > 6
            [minInd,~] = EOPnPL_Pm_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn); %返回最小重投影误差及其索引
    %         R0 = Rot(:, :, minInd);
    %         t0 = trans(:, minInd); 
            s0 = cgrs(:,minInd);
        else
            [minInd,~] = EOPnPL_Pm_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn); %返回最小重投影误差及其索引
            s0 = cgrs(:,minInd);
        end
        
% damped newton
        [s] = EOPnPL_Pm_damped_newton(s0,K,W,V,Xw,Ximgn,xs, xe, Xs, Xe);

        % compute the final R and t
        for i = 1 : size(s,2)
            [R(:,:,i),t(:,i),~,~] = EOPnPL_Pm_s2Rt(s(:,i),W,V);
        end

end

