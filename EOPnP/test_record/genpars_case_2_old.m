syms a1 a2 real;
syms v11 v12 v13 v14 v15 v16 v17 v18 v19 real;
syms v21 v22 v23 v24 v25 v26 v27 v28 v29 real;
v1 = [v11;v12;v13;v14;v15;v16;v17;v18;v19];
v2 = [v21;v22;v23;v24;v25;v26;v27;v28;v29];

%分解奇异向量
v1_3 = v1(1:3,:);    
v1_6 = v1(4:6,:);
v1_9 = v1(7:9,:);    %将奇异向量每三行分开，每个部分分别对应r1,r2,r3

v2_3 = v2(1:3,:);
v2_6 = v2(4:6,:);
v2_9 = v2(7:9,:);

%构建三个r的等式，并利用约束
r1 = a1*v1_3 + a2*v2_3;
r2 = a1*v1_6 + a2*v2_6;
r3 = a1*v1_9 + a2*v2_9;
f1 = dot(r1,r1) -1;
f2 = dot(r2,r2) -1;
f3 = dot(r3,r3) -1;
f4 = dot(r1,r2) -0;
f5 = dot(r1,r3) -0;
f6 = dot(r2,r3) -0;

%求代价函数C_a1a2
C_a1a2 = (f1)^2 + (f2)^2 + (f3)^2 + (f4)^2 + (f5)^2 + (f6)^2;

%分别对a1，a2求偏导
g1 = diff(C_a1a2,a1);
g2 = diff(C_a1a2,a2);

%使用coeffs函数提取系数，这里按照论文的方法将a2看作常数
coeff1 = coeffs(g1,a1);
coeff2 = coeffs(g2,a1);

coeff11 = coeff1(1,4);
coeff12 = coeff1(1,3);
coeff13 = coeff1(1,2);
coeff14 = coeff1(1,1);

coeff21 = coeff2(1,4);
coeff22 = coeff2(1,3);
coeff23 = coeff2(1,2);
coeff24 = coeff2(1,1);

%构造文中（22）矩阵
syms c11 c12 c13 c14 c21 c22 c23 c24
mat1 = [c11   0   0 c21   0   0;
        c12 c11   0 c22 c21   0;
        c13 c12 c11 c23 c22 c21;
        c14 c13 c12 c24 c23 c22;
          0 c14 c13   0 c24 c23;
          0   0 c14   0   0 c24];

f_a = det(mat1);
f_a = subs(f_a,c11,coeff11);
f_a = subs(f_a,c12,coeff12);
f_a = subs(f_a,c13,coeff13);
f_a = subs(f_a,c14,coeff14);
f_a = subs(f_a,c21,coeff21);
f_a = subs(f_a,c22,coeff22);
f_a = subs(f_a,c23,coeff23);
f_a = subs(f_a,c24,coeff24);
% matlabFunction(f_a)
[cfa,sfa] = coeffs(f_a,a2);
[cg1,sg1] = coeffs(g1,a1);








