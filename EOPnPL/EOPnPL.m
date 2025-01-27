function [R,t] = EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe)

% %normalize lines
% Xd = Xe - Xs;
% if (size(Xd, 2) > 0)
%     Xdn = sqrt(Xd(1, :).^2 + Xd(2, :).^2 + Xd(3, :).^2);
%     Xd = Xd ./ repmat(Xdn, 3, 1);
%     Xe = Xs + Xd;
% end

% 仅有线约束
if (isempty(Xw) && isempty(Ximgn))
    nl = size(Xs,2);
     if( all(abs(Xs(3,:) - repmat(Xe(3,1),1,nl)) < 1e-5) && ...
        all(abs(Xe(3,:) - repmat(Xe(3,1),1,nl)) < 1e-5)    ) 

        [R,t] = Planar_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);

     else

         [R,t] = EOPnPL_ordinary(Xw,Ximgn,xs,xe,Xs,Xe);

     end

%仅有点约束
elseif (isempty(Xs) && isempty(Xe))
    np = size(Xw,2);
     if( all(abs(Xw(3,:) - repmat(Xw(3,1),1,np)) < 1e-5) ) 
        
        [R,t] = Planar_EOPnPL(Xw,Ximgn,xs,xe,Xs,Xe);
     else
         tic;
         [R,t] = EOPnPL_ordinary(Xw,Ximgn,xs,xe,Xs,Xe);
         toc;
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
        [R,t] = EOPnPL_ordinary(Xw,Ximgn,xs,xe,Xs,Xe);
    end
end
