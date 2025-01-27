function [R,t] = EOPnPL_BJ_ordinary_1(Xw,Ximgn,xs,xe,Xs,Xe,iter1,iter2,PnP_iter)

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

    % 此处开始在模型直线上寻找与图像端点对应的点并使用 EOPnP 框架求解
    for i = 2:PnP_iter+1       
        [P1w_hat_temp] = EOPnPL_BJ_find_P1w_hat(R_temp,t_temp,xs,xe,Xs,Xe,nlines);       
        % EOPnP 框架或 EOPnL 框架
        [Rot,trans,cgrs,K,W,V] = EOPnPL_BJ_EOPnP_main([Xw P1w_hat_temp],[Ximgn [xs,xe]]);
%         [Rot, trans,cgrs,K,W,V] = EOPnPL_BJ_main(Xw,Ximgn,xs,xe,P1w_hat_temp(:,1:nlines),P1w_hat_temp(:,nlines+1:2*nlines)); 
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






