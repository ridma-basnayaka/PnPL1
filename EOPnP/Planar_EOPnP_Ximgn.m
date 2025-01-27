function [R,t] = Planar_EOPnP_Ximgn(Xw,Ximgn)
%输入参数要求：
% Ximgn ----- 归一化图像点坐标矩阵，每列为一个点的坐标，坐标要求非齐次（2*n）
% Xw ----- 点世界坐标矩阵，每列为一个点的坐标，坐标要求非齐次（3*n）
    n = size(Xw,2);
    A = [800 0 0;
         0 800 0;
         0 0 1];
     Ximg = A * [Ximgn;ones(1,n)];
     Ximg = Ximg(1:2,:);
     
    [K,W,V] = Planar_compute_K_Ximgn(Xw,Ximgn);

    [~,~,Vk] = svd(K,0);

    v1 = Vk(:,6);
    v2 = Vk(:,5);
    r = [Planar_case_1(v1) Planar_case_2_new(v1,v2)];

    %恢复完整的r
    c1 = r(1:2:5,:);    %第一列
    c2 = r(2:2:6,:);    %第二列
    c3 = cross(c1,c2);  %第三列
    rcomplete = [r(1:2,:);c3(1,:);r(3:4,:);c3(2,:);r(5:6,:);c3(3,:)];    

%% 恢复R和t
        [s_all,R_result_all,t_all] = Planar_recover_all_at_once(rcomplete,W,V);
        %除去相机平面后面的位姿结果
        ind = (t_all(3,:) > 0);
        s_all = s_all(:,ind);
        R_result_all = R_result_all(:,:,ind); 
        t_all = t_all(:,ind);


        for i = 1 : size(s_all,2)
            R = R_result_all(:,:,i);
            t = t_all(:,i);
%             cost_val(index) = Planar_cost_function_val(s,K);
            erri = reprojection_error_usingRT(Xw',Ximg',R,t,A);
            error(i) = erri; 
        end
        if n > 6
            [~,ind_min] = min(abs(error));
            s0 = s_all(:,ind_min);
        else
            s0 = s_all;
        end
        
        
    % damped newton
    s = Planar_damped_newton(s0,K,W,V,Xw,Ximg,A);
    
    % compute the final R and t
    for i = 1 : size(s,2)
            [R(:,:,i),t(:,i),~,~] = Planar_s2Rt(s(:,i),W,V);
    end


end


%     index = 1;
%     for i = 1 : size(rcomplete,2)
%         [R,t,s] = Planar_recover_try( rcomplete(:,i) , W , V);
%         proj = R*Xw + t*ones(1,n);
%         if isempty(R) || (t(3) < 0) || (min(proj(3,:)) < 0)
%             continue;
%         end
%         erri = reprojection_error_usingRT(Xw',Ximg',R,t,A);
%         error(index) = erri; ss(:,index) = s;
% %             Rot(:,:,index) = R; trans(:,index) = t; 
%         index = index + 1;
%     end
% 
%     if isempty(error)
%         R = []; t = []; err = [];
%         return;
%     end
%     [err,ind_min] = min(abs(error));
%     s0 = ss(:,ind_min);
