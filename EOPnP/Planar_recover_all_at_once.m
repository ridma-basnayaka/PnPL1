function [s_all,R_result_all,t_all] = Planar_recover_all_at_once(r,W,V)
n = size(r,2);
R_all = zeros(3,3,n);
for i = 1 : n
    R_all(:,:,i) = reshape(r(:,i),3,3)';
    [Ur,~,Vr] = svd( R_all(:,:,i) );
    R_all(:,:,i) = Ur*Vr';
end

%recover s
A = R_all(1,2,:) - R_all(2,1,:);
B = R_all(1,3,:) - R_all(3,1,:);
C = R_all(2,3,:) - R_all(3,2,:);

%第一组可能的s
s11 = -(2.*C + C.*(- A.^2 - B.^2 - C.^2 + 4).^(1/2))./(A.^2 + B.^2 + C.^2);
s21 = (2.*B + B.*(- A.^2 - B.^2 - C.^2 + 4).^(1/2))./(A.^2 + B.^2 + C.^2);
s31 = -(2.*A + A.*(- A.^2 - B.^2 - C.^2 + 4).^(1/2))./(A.^2 + B.^2 + C.^2);
s_all_1 = [s11;s21;s31];

%第二组可能的s
s12 = -(2.*C - C.*(- A.^2 - B.^2 - C.^2 + 4).^(1/2))./(A.^2 + B.^2 + C.^2);
s22 = (2.*B - B.*(- A.^2 - B.^2 - C.^2 + 4).^(1/2))./(A.^2 + B.^2 + C.^2);
s32 = -(2.*A - A.*(- A.^2 - B.^2 - C.^2 + 4).^(1/2))./(A.^2 + B.^2 + C.^2);
s_all_2 = [s12;s22;s32];

f1_all = R_all(1, 1,:) + R_all(2, 2,:) + R_all(3, 3,:) - (3 - s11 .* s11 - s21 .* s21 - s31 .* s31) ./ (1 + s11 .* s11 + s21 .* s21 + s31 .* s31);
f2_all = R_all(1, 1,:) + R_all(2, 2,:) + R_all(3, 3,:) - (3 - s12 .* s12 - s22 .* s22 - s32 .* s32) ./ (1 + s12 .* s12 + s22 .* s22 + s32 .* s32);

ind = (abs(f1_all) < abs(f2_all));
s_all = s_all_1.*(ind) + s_all_2.*(~ind);

% s2Rt
s1 = s_all(1,1,:);    s2 = s_all(2,1,:);    s3 = s_all(3,1,:);
Rbar_all = ...
   [s1.*s1 - s2.*s2 - s3.*s3 + 1,  2.*s1.*s2 - 2.*s3, 2.*s2 + 2.*s1.*s3;
            2.*s3 + 2.*s1.*s2, - s1.*s1 + s2.*s2 - s3.*s3 + 1, 2.*s2.*s3 - 2.*s1;
            2.*s1.*s3 - 2.*s2,  2.*s1 + 2.*s2.*s3, - s1.*s1 - s2.*s2 + s3.*s3 + 1];
coeff = (1+s1.*s1+s2.*s2+s3.*s3);   
R_result_all = Rbar_all./ coeff;

rbar_all = reshape( permute(Rbar_all,[2,1,3]) , 9 ,n);  % permute函数用于转置前两维

T_all = pinv(W) * V * rbar_all([1 2 4 5 7 8],:);
coeff = reshape(coeff,1,n);      %转为二维数组
t_all = T_all ./ coeff;
s_all = reshape(s_all,3,n);




end

