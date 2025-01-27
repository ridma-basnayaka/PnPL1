function [R2,t2] = EOPnPL_BJ_ordinary_2(Xw,Ximgn,xs,xe,Xs,Xe,iter1,iter2,PnP_iter)

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

        % 此处开始在模型直线上寻找与图像端点对应的点并使用EOPnP框架求解
        for i = 2:PnP_iter+1           
            [P1w_hat_temp] = EOPnPL_BJ_find_P1w_hat(R_temp,t_temp,xs,xe,Xs,Xe,nlines);       
            [Rot,trans,cgrs,K,W,V] = EOPnPL_BJ_EOPnP_main([Xw P1w_hat_temp],[Ximgn [xs,xe]]);
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
        end

        % 使用 PnPL 里面的方法对模型进行更新
        Xsf = R_temp*Xs + repmat(t_temp, 1, nlines);
        Xef = R_temp*Xe + repmat(t_temp, 1, nlines);               
        LcUnNormed = cross(Xsf, Xef);
        LcNorms = sqrt(LcUnNormed(1,:).^2 + LcUnNormed(2,:).^2);
        Lcf = LcUnNormed ./ repmat(LcNorms, 3, 1);
        [Xsc, Xec] = moveLineCloser2DLcVector(Xsf,Xef, xs, xe, Lcf); % 输出 Xsl、Xel 是相机坐标系坐标，函数外应转换为世界坐标系坐标
        Xs = R_temp'*Xsc - repmat(R_temp'*t_temp, 1, nlines);
        Xe = R_temp'*Xec - repmat(R_temp'*t_temp, 1, nlines);
        % 此处可选择是否再求解一次，还是直接迭代优化
        % [Rot,trans,cgrs,K,W,V] = EOPnPL_BJ_main(Xw,Ximgn,xs,xe,Xs,Xe);
        [K,W,V] = EOPnPL_BJ_compute_K(Xw,Ximgn,xs,xe,Xs,Xe);
    end  


[minInd,~] = EOPnPL_BJ_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn); %返回最小重投影误差及其索引 
R1 = Rot(:, :, minInd);
t1 = trans(:, minInd); 

s = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe, 15);

[R2,t2] = s2Rt(s,W,V);


end






