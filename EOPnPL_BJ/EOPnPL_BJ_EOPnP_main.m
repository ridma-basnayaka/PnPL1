function [Rot,trans,cgrs,K,W,V] = EOPnPL_BJ_EOPnP_main(Xw,Ximgn)
    
        [K,W,V] = EOPnPL_BJ_compute_K(Xw,Ximgn,[],[],[],[]);        
        [~,~,Vk] = svd(K,0);
        v1 = Vk(:,9);  v2 = Vk(:,8);  v3 = Vk(:,7);  v4 = Vk(:,6);
        r = [EOPnPL_BJ_case_1(v1), EOPnPL_BJ_case_2_new(v1,v2), EOPnPL_BJ_case_3_new(v1,v2,v3), EOPnPL_BJ_case_4_new(v1,v2,v3,v4)];

        %  recover R and t   
        [cgrs,Rot,trans] = EOPnPL_BJ_recover_all_at_once(r,W,V);
        if isempty(Rot)
            Rot = []; trans = []; 
            return;
        end
        
%         % 1.先选出重投影误差最小的解，再进行迭代优化
%         [minInd,~] = EOPnPL_BJ_Reproj_err(Rot, trans, xs, xe, Xs, Xe, Xw, Ximgn);
%         
%         % 2.先对所有解迭代优化，再选出重投影误差最小的解
end
