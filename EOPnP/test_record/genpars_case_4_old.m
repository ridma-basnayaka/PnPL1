syms a1 a2 a3 a4 k1 k2 k3 real
syms v11 v12 v13 v14 v15 v16 v17 v18 v19 real;
syms v21 v22 v23 v24 v25 v26 v27 v28 v29 real;
syms v31 v32 v33 v34 v35 v36 v37 v38 v39 real;
syms v41 v42 v43 v44 v45 v46 v47 v48 v49 real;
v1 = [v11;v12;v13;v14;v15;v16;v17;v18;v19];
v2 = [v21;v22;v23;v24;v25;v26;v27;v28;v29];
v3 = [v31;v32;v33;v34;v35;v36;v37;v38;v39];
v4 = [v41;v42;v43;v44;v45;v46;v47;v48;v49];

%分解奇异向量
v1_3 = v1(1:3,:);    
v1_6 = v1(4:6,:);
v1_9 = v1(7:9,:);    %将奇异向量每三行分开，每个部分分别对应r1,r2,r3

v2_3 = v2(1:3,:);
v2_6 = v2(4:6,:);
v2_9 = v2(7:9,:);

v3_3 = v3(1:3,:);
v3_6 = v3(4:6,:);
v3_9 = v3(7:9,:);

v4_3 = v4(1:3,:);
v4_6 = v4(4:6,:);
v4_9 = v4(7:9,:);

%旋转矩阵约束，顺便消去a1,a2,a3,a4
r1 = v1_3 + k1*v2_3 + k2*v3_3 + k3*v4_3;
r2 = v1_6 + k1*v2_6 + k2*v3_6 + k3*v4_6;
r3 = v1_9 + k1*v2_9 + k2*v3_9 + k3*v4_9;

%先考虑后三个约束，得到矩阵C元素
f4 = dot(r1,r2);
f5 = dot(r1,r3);
f6 = dot(r2,r3);

%提取矩阵M和N_k的系数
f = [f4;f5;f6];

f4_k3 = coeffs(f(1,1),k3);
f4_k3k2 = coeffs(f4_k3(1,1),k2);
f4_k2 = coeffs(f(1,1),k2);
f4_k2k3 = coeffs(f4_k2(1,1),k3);
f4_k2k3_2 = coeffs(f4_k2(1,2),k3);

f5_k3 = coeffs(f(2,1),k3);
f5_k3k2 = coeffs(f5_k3(1,1),k2);
f5_k2 = coeffs(f(2,1),k2);
f5_k2k3 = coeffs(f5_k2(1,1),k3);
f5_k2k3_2 = coeffs(f5_k2(1,2),k3);

f6_k3 = coeffs(f(3,1),k3);
f6_k3k2 = coeffs(f6_k3(1,1),k2);
f6_k2 = coeffs(f(3,1),k2);
f6_k2k3 = coeffs(f6_k2(1,1),k3);
f6_k2k3_2 = coeffs(f6_k2(1,2),k3);

M = [f4_k3k2(1,3) f4_k2k3(1,3) f4_k2k3_2(1,2);
     f5_k3k2(1,3) f5_k2k3(1,3) f5_k2k3_2(1,2);
     f6_k3k2(1,3) f6_k2k3(1,3) f6_k2k3_2(1,2)];
 
N_k = (-1) * [f4_k3k2(1,2) f5_k3k2(1,2) f6_k3k2(1,2);
              f4_k2k3(1,2) f5_k2k3(1,2) f6_k2k3(1,2);
              f4_k3k2(1,1) f5_k3k2(1,1) f6_k3k2(1,1)]';
          
Tk = pinv(M) * N_k;











