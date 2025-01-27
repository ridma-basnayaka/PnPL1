function [Xs_new,Xe_new] =EOPnPL_BJ_adjust_length_new_projected_model_line(xe,xs,P1c_hat,R0,t0,nlines)
% 对上述方法的进一步修改：使新的模型直线投影的长度与图像直线长度接近

length_image_lines = sum((xe-xs).^2);
Xs_new_c = P1c_hat(:,1:nlines);  Xe_new_c = P1c_hat(:,nlines+1:2*nlines);
Xs_new_projected = Xs_new_c(1:2,:) ./ Xs_new_c(3,:);
Xe_new_projected = Xe_new_c(1:2,:) ./ Xe_new_c(3,:);
length_projected_model_lines = sum((Xe_new_projected-Xs_new_projected).^2);
d1 = length_projected_model_lines - length_image_lines; % 图像直线长度和新模型直线投影的长度之差
% 长度长了往里缩，短了往外拉，两点各动二分之一的距离
V_new_projected_Xs_Xe = Xe_new_projected - Xs_new_projected;
Xs_new_projected = Xs_new_projected + 0.5*d1.*(V_new_projected_Xs_Xe);
Xe_new_projected = Xe_new_projected - 0.5*d1.*(V_new_projected_Xs_Xe);
% 将微调后的模型直线投影端点再重新反投影回相机坐标系和世界坐标系，用的是文档（2）的式（4）
Vc_xy = repmat(Xe_new_c(1:2,:) - Xs_new_c(1:2,:),1,2);
det_K = Vc_xy(1,:).*[Xs_new_projected(2,:),Xe_new_projected(2,:)]...
            - Vc_xy(2,:).*[Xs_new_projected(1,:),Xe_new_projected(1,:)];
nominator = Vc_xy(1,:).*P1c_hat(2,:) - Vc_xy(2,:).*P1c_hat(1,:);
depth = nominator ./ det_K;
P1c_hat = [Xs_new_projected,Xe_new_projected;ones(1,2*nlines)] .* depth;
P1w_hat_temp = R0\(P1c_hat - repmat(t0,1,2*nlines));
Xs_new = P1w_hat_temp(:,1:nlines);  Xe_new = P1w_hat_temp(:,nlines+1:2*nlines);

end

