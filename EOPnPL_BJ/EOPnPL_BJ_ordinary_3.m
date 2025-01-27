function [R,t] = EOPnPL_BJ_ordinary_3(Xw,Ximgn,xs,xe,Xs,Xe,iter1,iter2,PnP_iter)

nlines = size(Xs,2);
[Rot, trans,cgrs,K,W,V] = EOPnPL_BJ_main(Xw,Ximgn,xs,xe,Xs,Xe);

if nlines > 0
    R = zeros(3,3,PnP_iter+1);  t = zeros(3,PnP_iter+1);  
    [minInd,~] = EOPnPL_BJ_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
    
    if iter1 == 0
        R_temp = Rot(:, :, minInd);  t_temp = trans(:, minInd);
        R(:,:,1) = R_temp;  t(:,1) = t_temp;                
    else
        s_temp = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe,iter1);
        [R_temp,t_temp] = EOPnPL_BJ_s2Rt(s_temp,W,V);
        R(:,:,1) = R_temp;  t(:,1) = t_temp;
    end

    % 此处开始使用 PnPL 里面的方法对模型进行更新并使用EOPnP框架求解
    for i = 2:PnP_iter+1 
        % PnPL 里更新模型的方法
        Xsf = R_temp*Xs + repmat(t_temp, 1, nlines);
        Xef = R_temp*Xe + repmat(t_temp, 1, nlines);               
        LcUnNormed = cross(Xsf, Xef);
        LcNorms = sqrt(LcUnNormed(1,:).^2 + LcUnNormed(2,:).^2);
        Lcf = LcUnNormed ./ repmat(LcNorms, 3, 1);
        [Xsc, Xec] = moveLineCloser2DLcVector(Xsf,Xef, xs, xe, Lcf); % 输出 Xsl、Xel 是相机坐标系坐标，函数外应转换为世界坐标系坐标
        Xs = R_temp'*Xsc - repmat(R_temp'*t_temp, 1, nlines);
        Xe = R_temp'*Xec - repmat(R_temp'*t_temp, 1, nlines);
        
        % 使用 EOPnP 或 EOPnPL 框架求解
        [Rot,trans,cgrs,K,W,V] = EOPnPL_BJ_EOPnP_main([Xw [Xs,Xe]],[Ximgn [xs,xe]]);
%         [Rot, trans,cgrs,K,W,V] = EOPnPL_BJ_main(Xw,Ximgn,xs,xe,Xs,Xe);
        [minInd,~] = EOPnPL_BJ_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
        
        % 1.只对最后一次循环的结果进行优化  
        if (iter2 ~= 0) && (i == PnP_iter+1)  
            s_temp = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe,iter2);
            [R_temp,t_temp] = EOPnPL_BJ_s2Rt(s_temp,W,V);
            R(:,:,i) = R_temp;  t(:,i) = t_temp;
        else  % iter2 == 0 或者未到最后一次循环
            R_temp = Rot(:, :, minInd);  t_temp = trans(:, minInd);
            R(:,:,i) = R_temp;  t(:,i) = t_temp;            
        end
        
%         % 2.每次循环的结果都优化
%         if iter2 == 0
%             R_temp = Rot(:, :, minInd);  t_temp = trans(:, minInd);
%             R(:,:,i) = R_temp;  t(:,i) = t_temp;
%         else  
%             s_temp = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe,iter2);
%             [R_temp,t_temp] = EOPnPL_BJ_s2Rt(s_temp,W,V);
%             R(:,:,i) = R_temp;  t(:,i) = t_temp;
%         end
    end
end  

R = R(:,:,end);  t = t(:,end);

end






