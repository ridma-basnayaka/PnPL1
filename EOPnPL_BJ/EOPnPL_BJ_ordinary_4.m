function [R,t] = EOPnPL_BJ_ordinary_4(Xw,Ximgn,xs,xe,Xs,Xe,iter1,iter2,PnP_iter)


nlines = size(Xs,2);
[Rot, trans,cgrs,K,W,V] = EOPnPL_BJ_main(Xw,Ximgn,xs,xe,Xs,Xe);

if nlines > 0
    [minInd,~] = EOPnPL_BJ_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
    if iter1 == 0
        R_temp = Rot(:, :, minInd);  t_temp = trans(:, minInd);
    else
        s_temp = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe,iter1);
        [R_temp,t_temp] = EOPnPL_BJ_s2Rt(s_temp,W,V);
    end
    [P1w_hat_temp,P1c_hat_temp] = EOPnPL_BJ_find_P1w_hat(R_temp,t_temp,xs,xe,Xs,Xe,nlines); 
    % Xs_new = P1w_hat_temp(:,1:nlines);  Xe_new = P1w_hat_temp(:,nlines+1:2*nlines);

    % 对上述方法的进一步修改：使新的模型直线投影的长度与图像直线长度接近
    %（没效果，上一步产生的新模型直线的投影的长度本就与图像直线长度很接近）
    [Xs_new,Xe_new] =EOPnPL_BJ_adjust_length_new_projected_model_line(xe,xs,P1c_hat_temp,R_temp,t_temp,nlines);

    [K,W,V] = EOPnPL_BJ_compute_K(Xw,Ximgn,xs,xe,Xs_new,Xe_new);

end

s = EOPnPL_BJ_damped_newton(cgrs(:,minInd),K,W,V,Xw,Ximgn,xs, xe, Xs, Xe, 5);
[R,t] = s2Rt(s,W,V);


end






