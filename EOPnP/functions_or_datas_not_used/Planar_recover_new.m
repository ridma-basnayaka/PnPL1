function [R,t,s] = Planar_recover_new(r,W,V)
%using t calculated before to recover 
%rotation R and translation t

R = reshape(r,3,3)';
[Ur,~,Vr] = svd(R);
R_e = Ur*Vr';

%recover s
A = R_e(1,2) - R_e(2,1);
B = R_e(1,3) - R_e(3,1);
C = R_e(2,3) - R_e(3,2);
%第一组可能的s
s11 = -(2*C + C*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
s21 = (2*B + B*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
s31 = -(2*A + A*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
f = R_e(1, 1) + R_e(2, 2) + R_e(3, 3) - (3 - s11 * s11 - s21 * s21 - s31 * s31) / (1 + s11 * s11 + s21 * s21 + s31 * s31);
if (abs(f) < 1e-5)
    s = [s11;s21;s31];
else
    %第二组s
    s12 = -(2*C - C*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
    s22 = (2*B - B*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
    s32 = -(2*A - A*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
    f = R_e(1, 1) + R_e(2, 2) + R_e(3, 3) - (3 - s12 * s12 - s22 * s22 - s32 * s32) / (1 + s12 * s12 + s22 * s22 + s32 * s32);
    if (abs(f) < 1e-5)
        s = [s12;s22;s32];
    else
        s = []; %两组s都不符合，s无解置空。
    end
end



%recover t
if isempty(s)
    R = [];
    t = [];
    return;
else
    [R,t,~,~] = Planar_s2Rt(s,W,V);
end
   
end
