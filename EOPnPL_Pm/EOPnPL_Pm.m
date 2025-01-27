function [R,t] = EOPnPL_Pm(Xw,Ximgn,xs,xe,Xs,Xe)

if (isempty(Xw) && isempty(Ximgn))
    nl = size(Xs,2);
     if( all(abs(Xs(3,:) - repmat(Xe(3,1),1,nl)) < 1e-5) && ...
        all(abs(Xe(3,:) - repmat(Xe(3,1),1,nl)) < 1e-5)    ) 

        [R,t] = Planar_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);

     else

         [R,t] = EOPnPL_Pm_ordinary(Xw,Ximgn,xs,xe,Xs,Xe);

     end

%仅有点约束
elseif (isempty(Xs) && isempty(Xe))
    np = size(Xw,2);
     if( all(abs(Xw(3,:) - repmat(Xw(3,1),1,np)) < 1e-5) ) 
        
        [R,t] = Planar_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);
     else
         [R,t] = EOPnPL_ordinary_Pm(Xw,Ximgn,xs,xe,Xs,Xe);
     end

%点线均有    
else
    np = size(Xw,2); %参考点数量
    nl = size(Xs,2);
    if( all(abs(Xw(3,:) - repmat(Xw(3,1),1,np)) < 1e-5) && ...
        all(abs(Xs(3,:) - repmat(Xw(3,1),1,nl)) < 1e-5) && ...
        all(abs(Xe(3,:) - repmat(Xw(3,1),1,nl)) < 1e-5)    ) 
        
        [R,t] = Planar_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);
        
    else
% 普通配置    
        [R,t] = EOPnPL_ordinary_Pm(Xw,Ximgn,xs,xe,Xs,Xe);
    end

end
