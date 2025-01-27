function [R,t,s] = recover_try(r,W,V,K)
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

cost_val_1 = cost_function_val(s1,K);
cost_val_2 = cost_function_val(s2,K);

%recover t
if (cost_val_1 < cost_val_2)
    s = s1;
    [R,t,~,~] = s2Rt(s,W,V);
else
    s = s2;
    [R,t,~,~] = s2Rt(s,W,V);
end
   
end
