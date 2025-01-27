function [Rot,trans,cgrs,K,W,V] = EOPnPL_Pm_main(Xw,Ximgn,xs,xe,Xs,Xe)

        [K,W,V] = EOPnPL_Pm_compute_K(Xw,Ximgn,xs,xe,Xs,Xe);%非共面情况
        
        [~,~,Vk] = svd(K,0);

        v1 = Vk(:,9);
        v2 = Vk(:,8);
        v3 = Vk(:,7);
        v4 = Vk(:,6);
        r = [EOPnPL_Pm_case_1(v1), EOPnPL_Pm_case_2_new(v1,v2), ...
            EOPnPL_Pm_case_3_new(v1,v2,v3), EOPnPL_Pm_case_4_new(v1,v2,v3,v4)];


        %  recover R and t   
        Rot = []; trans = [];
        [cgrs,Rot,trans] = EOPnPL_Pm_recover_all_at_once(r,W,V);

        if isempty(Rot)
            Rot = []; trans = []; 
            return;
        end

end



