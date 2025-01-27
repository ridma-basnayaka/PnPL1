function [R,t] = EOPnP(Xw,Ximgn)
%输入参数要求：
% Ximgn ----- 归一化图像点坐标矩阵，每列为一个点的坐标，坐标要求非齐次（2*n）
% Xw ----- 点世界坐标矩阵，每列为一个点的坐标，坐标要求非齐次（3*n）

     n = size(Xw,2); %参考点数量
%      A = [800 0 0; 0 800 0; 0 0 1];
%      Ximg = A * [Ximgn;ones(1,n)];
%      Ximg = Ximg(1:2,:);

     
    if( all(abs(Xw(3,:) - repmat(Xw(3,1),1,n)) < 1e-5) ) %共面情况
        
        [R,t] = Planar_EOPnP_Ximgn(Xw,Ximgn);        
    else
% tic;
        [K,W,V] = compute_K_Ximgn(Xw,Ximgn);%非共面情况
% toc;
% tic;
        [~,~,Vk] = svd(K,0);
% toc;

        v1 = Vk(:,9);
        v2 = Vk(:,8);
        v3 = Vk(:,7);
        v4 = Vk(:,6);

        r = [case_1(v1), case_2_new(v1,v2), case_3_new(v1,v2,v3), case_4_new(v1,v2,v3,v4)];


        [s_all,R_result_all,t_all] = recover_all_at_once(r,W,V);
        
        
        %不去除相机平面后面的位姿结果
        s_all = s_all;
        R_result_all = R_result_all; 
        t_all = t_all;        

        for i = 1 : size(s_all,2)
            R = R_result_all(:,:,i);
            t = t_all(:,i);
%             cost_val(index) = cost_function_val(s,K);
            erri = reprojection_error_usingRT(Xw',Ximgn',R,t,eye(3,3)); % 直接使用归一化坐标计算重投影误差
            error(i) = erri; 
        end
        if n > 6
            [~,ind_min] = min(abs(error));
            s0 = s_all(:,ind_min);
        else
            s0 = s_all;
        end    

        % damped newton
        s = damped_newton(s0,K,W,V,Xw,Ximgn,eye(3,3));

        % compute the final R and t
        for i = 1 : size(s,2)
            [R(:,:,i),t(:,i),~,~] = s2Rt(s(:,i),W,V);
        end
     
    end

end



% version 2
%         index = 1;
%         %recover R and t
%         for i = 1 : size(r,2)
%             [R,t,s] = recover_try2( r(:,i) , W , V);
%             proj = R*Xw + t*ones(1,n);
%             if isempty(R) || (t(3) < 0) || (min(proj(3,:)) < 0)
%                 continue;
%             end
%             erri = reprojection_error_usingRT(Xw',Ximg',R,t,A);
%             error(index) = erri; ss(:,index) = s;
% %             Rot(:,:,index) = R; trans(:,index) = t; 
%             index = index + 1;
%         end
%         
%         if isempty(error)
%             R = []; t = []; err = [];
%             return;
%         end
%         [err,ind_min] = min(abs(error));
%         s0 = ss(:,ind_min);


% version 1
% n1 = size(r,2);
%         err = inf;
%         for i = 1 : n1
%             [Ri,ti,si] = recover_try2( r(:,i) , W , V);
% %             [Ri,ti,si] = recover_try( r(:,i) , W , V, K );
% %            [Ri,ti,si] = recover_new( r(:,i) , W , V );
%             if isempty(Ri) || (ti(3) < 0)
%                 continue;
%             else   
% %                 erri = reprojection_error_usingRT(Xw',Ximg',Ri,ti,A);
%                 erri = cost_function_val(si,K);
% %                 i
% %                 erri
%                 if (erri < err)
%                     err = erri;
%                     R0 = Ri;
%                     t0 = ti;
%                     s0 = si;
% %                     i
%                 end
%             end
%         end
% %         rn = reshape(R0',9,1);
% %         [R0n,t0n,s0n] = recover_try2( -rn , W , V ); 
% %         cost_val_1 = cost_function_val(s0,K);
% %         cost_val_2 = cost_function_val(s0n,K);
% %         if(cost_val_1 > cost_val_2)
% %             s0 = s0n;
% %         end
