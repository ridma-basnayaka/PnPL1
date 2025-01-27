function [R,t,s] = recover_try2(r,W,V)
%using t calculated before to recover 
%rotation R and translation t
if isempty(r)
    R = [];
    t = [];
    return;
end

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
s1 = [s11;s21;s31];

%第二组可能的s
s12 = -(2*C - C*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
s22 = (2*B - B*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
s32 = -(2*A - A*(- A^2 - B^2 - C^2 + 4)^(1/2))/(A^2 + B^2 + C^2);
s2 = [s12;s22;s32];

f1 = R_e(1, 1) + R_e(2, 2) + R_e(3, 3) - (3 - s11 * s11 - s21 * s21 - s31 * s31) / (1 + s11 * s11 + s21 * s21 + s31 * s31);
f2 = R_e(1, 1) + R_e(2, 2) + R_e(3, 3) - (3 - s12 * s12 - s22 * s22 - s32 * s32) / (1 + s12 * s12 + s22 * s22 + s32 * s32);

%recover t
if (abs(f1) < abs(f2))
    s = s1;
    [R,t,~,~] = s2Rt(s,W,V);
else
    s = s2;
    [R,t,~,~] = s2Rt(s,W,V);
end
   
end
