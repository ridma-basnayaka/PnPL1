function [Rot,trans,cgrs,K,W,V] = EOPnPL_main(Xw,Ximgn,xs,xe,Xs,Xe)

[K,W,V] = compute_K_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);%非共面情况

[~,~,Vk] = svd(K,0);

v1 = Vk(:,9);
v2 = Vk(:,8);
v3 = Vk(:,7);
v4 = Vk(:,6);
r = [case_1(v1), case_2_new(v1,v2), case_3_new(v1,v2,v3), case_4_new(v1,v2,v3,v4)];

%  recover R and t   
[cgrs,Rot,trans] = recover_all_at_once(r,W,V);

if isempty(Rot)
    Rot = []; trans = []; 
    return;
end

end




%         index = 1;
%         for i = 1 : size(r,2)
%             R = R_result_all(:,:,i);
%             t = t_all(:,i);
%             s = s_all(:,i);
%             if isempty(Xw)
%                 proj = R*[Xs Xe] + t*ones(1,2*nl);
%             elseif isempty(Xs)
%                  proj = R*[Xw] + t*ones(1,np);
%             else
%                 proj = R*[Xw Xs Xe] + t*ones(1,np + 2*nl);
%             end
%             if (t(3) < 0)  || (min(proj(3,:)) < 0)
%                 continue;
%             end
%             cgrs(:,index) = s;
%             Rot(:,:,index) = R; trans(:,index) = t;
%             index = index + 1;
%         end

