addpath others

%recover s
A = R(1,2) - R(2,1);
B = R(1,3) - R(3,1);
C = R(2,3) - R(3,2);

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

f1 = R(1, 1) + R(2, 2) + R(3, 3) - (3 - s11 * s11 - s21 * s21 - s31 * s31) / (1 + s11 * s11 + s21 * s21 + s31 * s31);
f2 = R(1, 1) + R(2, 2) + R(3, 3) - (3 - s12 * s12 - s22 * s22 - s32 * s32) / (1 + s12 * s12 + s22 * s22 + s32 * s32);

%recover t
if (abs(f1) < abs(f2))
    s = s1;
else
    s = s2;
end


disp(s);
disp(t);

[R1,t1] = EOPnPL(XXw,xxn,xs,xe,Xs,Xe);


index_best = 1;
error = inf;
y = 1e10*ones(2,1);
for jjj = 1:size(R1,3)
    tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
    if sum(tempy) < error
        y = tempy;
        error = sum(tempy);
        index_best = jjj;
    end
end
disp(y);