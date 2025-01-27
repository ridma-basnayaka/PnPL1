function [minInd,minErr] = Reproj_err(R0, t0, xs, xe, Xs, Xe, U, u)
% 输入R0为三维矩阵，存储待选的旋转矩阵R，第三维度数量为旋转矩阵数量
% 输入t0为二维矩阵，存储待选的平移向量t，第二维度数量为平移向量数量

    if (size(R0, 3) == 0)
        minInd = 0;
        return;        
    end
    minErr = 1e10;
    minInd = 1;
    Xw = [Xs, Xe, U];
    Xw_h=[Xw;ones(1,size(Xw,2))];    
    for i = 1:size(R0, 3)
        R = R0(:, :, i);
        t = t0(:, i); 
        
        Xc = [R,t] * Xw_h;
        if min(Xc(3,:)) < 0
            continue;
        end  %排除点在相机后面的结果
            
        if (size(U, 2) > 0)
            projs = R*U + repmat(t, 1, size(U, 2));
            projs(1:2, :) = projs(1:2, :) ./ repmat(projs(3, :), 2, 1);
            res = u(1:2, :) - projs(1:2, :);
            aggPtErr = mean(sqrt(res(1, :).^2 + res(2, :).^2));
        else
            aggPtErr = 0;
        end
               
        aggErrLi = 0;
        if (size(Xs, 2) > 0)
            Xsp = R*Xs + repmat(t, 1, size(Xs, 2));
            Xep = R*Xe + repmat(t, 1, size(Xe, 2));        

            LcUnNormed = cross(Xsp, Xep);
            LcNorms = sqrt(LcUnNormed(1,:).^2 + LcUnNormed(2,:).^2);
            LcNormed = LcUnNormed ./ repmat(LcNorms, 3, 1);

            aggErrLi = sum(abs(LcNormed(1,:).*xs(1,:) + LcNormed(2,:).*xs(2,:) + LcNormed(3,:)));
            aggErrLi = aggErrLi + sum(abs(LcNormed(1,:).*xe(1,:) + LcNormed(2,:).*xe(2,:) + LcNormed(3,:)));
        end
        
        errVal = aggPtErr + aggErrLi;
        if (errVal < minErr)
            minErr = errVal;
            minInd = i;
        end
    end    
end
