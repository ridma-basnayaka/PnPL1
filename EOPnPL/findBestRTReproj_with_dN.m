function [minInd,Xsf,Xef,Lcf,s] = findBestRTReproj_with_dN(R0, t0, xs, xe, Xs, Xe, U, u, cgrs, K, W, V)
% 输入R0为三维矩阵，存储待选的旋转矩阵R，第三维度数量为旋转矩阵数量
% 输入t0为二维矩阵，存储待选的平移向量t，第二维度数量为平移向量数量
% U为参考点世界坐标，u为图像点归一化像素坐标
% 输出Xsf、Xef、Lcf分别是选定外参下相机坐标系中对应直线、起始点、终止点的坐标
    if (size(R0, 3) == 0)
        minInd = 0;
        Xsf = [];
        Xef = [];
        Lcf = [];
        return;        
    end
    minInd = 1;
    errVal = [];
    for i = 1:size(R0, 3)
        R = R0(:, :, i);
        t = t0(:, i); 
        
        if (size(U, 2) > 0)
            projs = R*U + repmat(t, 1, size(U, 2));
            projs(1:2, :) = projs(1:2, :) ./ repmat(projs(3, :), 2, 1);
            res = u(1:2, :) - projs(1:2, :);
            aggPtErr = mean(sqrt(res(1, :).^2 + res(2, :).^2));
        else
            aggPtErr = 0;
        end
        
        Xsp = R*Xs + repmat(t, 1, size(Xs, 2));
        Xep = R*Xe + repmat(t, 1, size(Xe, 2));        
        
        LcUnNormed = cross(Xsp, Xep);
        LcNorms = sqrt(LcUnNormed(1,:).^2 + LcUnNormed(2,:).^2);
        LcNormed = LcUnNormed ./ repmat(LcNorms, 3, 1);
              
        aggErrLi = sum(abs(LcNormed(1,:).*xs(1,:) + LcNormed(2,:).*xs(2,:) + LcNormed(3,:)));
        aggErrLi = aggErrLi + sum(abs(LcNormed(1,:).*xe(1,:) + LcNormed(2,:).*xe(2,:) + LcNormed(3,:)));
        errVal(i) = aggPtErr + aggErrLi;
    end    
    [~,minInd] = sort(errVal);
    % damped newton
    s0 = cgrs(:,[minInd(1:5)]);
    s = Planar_damped_newton_EOPnPL(s0,K,W,V,U,u,xs, xe, Xs, Xe);
    for i = 1 : size(s,2)
        [R(:,:,i),t(:,i),~,~] = Planar_s2Rt(s(:,i),W,V);
    end
    [minInd2,~] = Reproj_err(R, t, xs, xe, Xs, Xe, U, u);
    R = R(:,:,minInd2);  t = t(:,minInd2);  
    minInd = minInd(minInd2);  s = s(:,minInd2);
    
    Xsp = R*Xs + repmat(t, 1, size(Xs, 2));
    Xep = R*Xe + repmat(t, 1, size(Xe, 2));        
    LcUnNormed = cross(Xsp, Xep);
    LcNorms = sqrt(LcUnNormed(1,:).^2 + LcUnNormed(2,:).^2);
    LcNormed = LcUnNormed ./ repmat(LcNorms, 3, 1);    
    Xsf = Xsp;
    Xef = Xep;
    Lcf = LcNormed;
end
