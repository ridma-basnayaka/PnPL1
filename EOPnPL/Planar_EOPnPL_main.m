function [Rot,trans,cgrs,K,W,V] = Planar_EOPnPL_main(Xw,Ximgn,xs,xe,Xs,Xe)

    np = size(Xw,2); %参考点数量  
    nl = size(Xs,2);
    [K,W,V] = Planar_compute_K_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);

    [~,~,Vk] = svd(K,0);

    v1 = Vk(:,6);
    v2 = Vk(:,5);
%     r = [Planar_case_1(v1) Planar_case_2_new(v1,v2)];
    r = [Planar_case_1(v1) Planar_case_2_new(v1,v2) Planar_case_3_new(v1,v2,Vk(:,4))];   
%     r = [Planar_case_1(v1) Planar_case_2_new(v1,v2) Planar_case_1_z(v1) Planar_case_2_z(v1,v2)];    
    
    

    %恢复完整的r
    c1 = r(1:2:5,:);    %第一列
    c2 = r(2:2:6,:);    %第二列
    c3 = cross(c1,c2);  %第三列
    rcomplete = [r(1:2,:);c3(1,:);r(3:4,:);c3(2,:);r(5:6,:);c3(3,:)];    

    % 恢复R和t
    Rot = []; trans = [];
    [s_all,R_result_all,t_all] = Planar_recover_all_at_once(rcomplete,W,V);
    ind = (t_all(3,:) > 0);
    cgrs = s_all(:,ind);
    Rot = R_result_all(:,:,ind); 
    trans = t_all(:,ind);

    if isempty(Rot)
        Rot = []; trans = []; 
        return;
    end

end


%         index = 1;
%         for i = 1 : size(rcomplete,2)
%             R = R_result_all(:,:,i);
%             t = t_all(:,i);
%             s = s_all(:,i);
% %             if isempty(Xw)
% %                 proj = R*[Xs Xe] + t*ones(1,2*nl);
% %             elseif isempty(Xs)
% %                  proj = R*[Xw] + t*ones(1,np);
% %             else
% %                 proj = R*[Xw Xs Xe] + t*ones(1,np + 2*nl);
% %             end
%             if (t(3) < 0) %|| (min(proj(3,:)) < 0)
%                 continue;
%             end
%             cgrs(:,index) = s;
%             Rot(:,:,index) = R; trans(:,index) = t; 
%             index = index + 1;
%         end
