function [] = case_2_twomethods( v1 , v2 )

%分解奇异向量
v11 = v1(1:3,:);    
v12 = v1(4:6,:);
v13 = v1(7:9,:);    %将奇异向量每三行分开，每个部分分别对应r1,r2,r3

v21 = v2(1:3,:);
v22 = v2(4:6,:);
v23 = v2(7:9,:);

%设置系数未知数
syms a1 a2 real;

%构建三个r的等式，并利用约束
r1 = a1.*v11 + a2.*v21;
r2 = a1.*v12 + a2.*v22;
r3 = a1.*v13 + a2.*v23;
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
c1 = coeffs(g1,a1);
c2 = coeffs(g2,a1);

c11 = c1(1,4);
c12 = c1(1,3);
c13 = c1(1,2);
c14 = c1(1,1);

c21 = c2(1,4);
c22 = c2(1,3);
c23 = c2(1,2);
c24 = c2(1,1);

%构造文中（22）矩阵
mat1 = [c11   0   0 c21   0   0;
        c12 c11   0 c22 c21   0;
        c13 c12 c11 c23 c22 c21;
        c14 c13 c12 c24 c23 c22;
          0 c14 c13   0 c24 c23;
          0   0 c14   0   0 c24];
      
%求上面矩阵的行列式，得到函数f_a
f_a = det(mat1);
[cfa_true,~] = coeffs(f_a,a2);
k9_true = vpa(cfa_true(1));
k7_true = vpa(cfa_true(2));
k5_true = vpa(cfa_true(3));
k3_true = vpa(cfa_true(4));
k1_true = vpa(cfa_true(5));


%其他方法
%一次性赋值用deal函数
% [v11;v12;v13;v14;v15;v16;v17;v18;v19] = deal(v1);
% [v21;v22;v23;v24;v25;v26;v27;v28;v29] = deal(v2);
%分解奇异向量
v1_3 = v1(1:3,:);    
v1_6 = v1(4:6,:);
v1_9 = v1(7:9,:);    %将奇异向量每三行分开，每个部分分别对应r1,r2,r3

v2_3 = v2(1:3,:);
v2_6 = v2(4:6,:);
v2_9 = v2(7:9,:);

A = [   dot(v1_3,v1_3)    2 * dot(v1_3,v2_3)                        dot(v2_3,v2_3);	
        dot(v1_6,v1_6)    2 * dot(v1_6,v2_6)                        dot(v2_6,v2_6);	
        dot(v1_9,v1_9)    2 * dot(v1_9,v2_9)                        dot(v2_9,v2_9);	
        dot(v1_3,v1_6)	( dot(v1_3,v2_6) + dot(v2_3,v1_6) )         dot(v2_3,v2_6);
        dot(v1_3,v1_9)	( dot(v1_3,v2_9) + dot(v2_3,v1_9) )         dot(v2_3,v2_9);
        dot(v1_6,v1_9)	( dot(v1_6,v2_9) + dot(v2_6,v1_9) )         dot(v2_6,v2_9) ];

A_new  = zeros(6,5);
for i = 1 : 6
    A_new(i,1) = A(i,1)^2;
    A_new(i,2) = 2*A(i,1)*A(i,2);
    A_new(i,3) = A(i,2)^2 + 2*A(i,1)*A(i,3);
    A_new(i,4) = 2*A(i,2)*A(i,3); 
    A_new(i,5) = A(i,3)^2;
end

Ans = ( A_new(1,:) + A_new(2,:) + A_new(3,:) + A_new(4,:) + A_new(5,:) + A_new(6,:) );
As = A(1,:) + A(2,:) + A(3,:);

k9 = 4096*Ans(1)^3*Ans(5)^3 - 3072*Ans(1)^2*Ans(2)*Ans(4)*Ans(5)^2 - 2048*Ans(1)^2*Ans(3)^2*Ans(5)^2 + 2304*Ans(1)^2*Ans(3)*Ans(4)^2*Ans(5) - 432*Ans(1)^2*Ans(4)^4 + 2304*Ans(1)*Ans(2)^2*Ans(3)*Ans(5)^2 - 96*Ans(1)*Ans(2)^2*Ans(4)^2*Ans(5) - 1280*Ans(1)*Ans(2)*Ans(3)^2*Ans(4)*Ans(5) + 288*Ans(1)*Ans(2)*Ans(3)*Ans(4)^3 + 256*Ans(1)*Ans(3)^4*Ans(5) - 64*Ans(1)*Ans(3)^3*Ans(4)^2 - 432*Ans(2)^4*Ans(5)^2 + 288*Ans(2)^3*Ans(3)*Ans(4)*Ans(5) - 64*Ans(2)^3*Ans(4)^3 - 64*Ans(2)^2*Ans(3)^3*Ans(5) + 16*Ans(2)^2*Ans(3)^2*Ans(4)^2;
k7 = - 12288*As(3)*Ans(1)^3*Ans(5)^2 + 6144*As(3)*Ans(1)^2*Ans(2)*Ans(4)*Ans(5) + 3072*As(2)*Ans(1)^2*Ans(2)*Ans(5)^2 + 4096*As(3)*Ans(1)^2*Ans(3)^2*Ans(5) - 2304*As(3)*Ans(1)^2*Ans(3)*Ans(4)^2 - 4608*As(2)*Ans(1)^2*Ans(3)*Ans(4)*Ans(5) + 4096*As(1)*Ans(1)^2*Ans(3)*Ans(5)^2 + 1728*As(2)*Ans(1)^2*Ans(4)^3 - 2304*As(1)*Ans(1)^2*Ans(4)^2*Ans(5) - 4608*As(3)*Ans(1)*Ans(2)^2*Ans(3)*Ans(5) + 96*As(3)*Ans(1)*Ans(2)^2*Ans(4)^2 + 192*As(2)*Ans(1)*Ans(2)^2*Ans(4)*Ans(5) - 2304*As(1)*Ans(1)*Ans(2)^2*Ans(5)^2 + 1280*As(3)*Ans(1)*Ans(2)*Ans(3)^2*Ans(4) + 1280*As(2)*Ans(1)*Ans(2)*Ans(3)^2*Ans(5) - 864*As(2)*Ans(1)*Ans(2)*Ans(3)*Ans(4)^2 + 2560*As(1)*Ans(1)*Ans(2)*Ans(3)*Ans(4)*Ans(5) - 288*As(1)*Ans(1)*Ans(2)*Ans(4)^3 - 256*As(3)*Ans(1)*Ans(3)^4 + 128*As(2)*Ans(1)*Ans(3)^3*Ans(4) - 1024*As(1)*Ans(1)*Ans(3)^3*Ans(5) + 192*As(1)*Ans(1)*Ans(3)^2*Ans(4)^2 + 864*As(3)*Ans(2)^4*Ans(5) - 288*As(3)*Ans(2)^3*Ans(3)*Ans(4) - 288*As(2)*Ans(2)^3*Ans(3)*Ans(5) + 192*As(2)*Ans(2)^3*Ans(4)^2 - 288*As(1)*Ans(2)^3*Ans(4)*Ans(5) + 64*As(3)*Ans(2)^2*Ans(3)^3 - 32*As(2)*Ans(2)^2*Ans(3)^2*Ans(4) + 192*As(1)*Ans(2)^2*Ans(3)^2*Ans(5) - 32*As(1)*Ans(2)^2*Ans(3)*Ans(4)^2;
k5 = 12288*Ans(5)*Ans(1)^3*As(3)^2 - 3072*Ans(1)^2*Ans(2)*Ans(4)*As(3)^2 - 6144*Ans(5)*Ans(1)^2*Ans(2)*As(2)*As(3) - 2048*Ans(1)^2*Ans(3)^2*As(3)^2 + 4608*Ans(1)^2*Ans(3)*Ans(4)*As(2)*As(3) - 8192*Ans(5)*Ans(1)^2*Ans(3)*As(1)*As(3) + 2048*Ans(5)*Ans(1)^2*Ans(3)*As(2)^2 + 2304*Ans(1)^2*Ans(4)^2*As(1)*As(3) - 2304*Ans(1)^2*Ans(4)^2*As(2)^2 + 3072*Ans(5)*Ans(1)^2*Ans(4)*As(1)*As(2) + 2304*Ans(1)*Ans(2)^2*Ans(3)*As(3)^2 - 192*Ans(1)*Ans(2)^2*Ans(4)*As(2)*As(3) + 4608*Ans(5)*Ans(1)*Ans(2)^2*As(1)*As(3) - 1280*Ans(1)*Ans(2)*Ans(3)^2*As(2)*As(3) - 2560*Ans(1)*Ans(2)*Ans(3)*Ans(4)*As(1)*As(3) + 640*Ans(1)*Ans(2)*Ans(3)*Ans(4)*As(2)^2 - 1536*Ans(5)*Ans(1)*Ans(2)*Ans(3)*As(1)*As(2) + 960*Ans(1)*Ans(2)*Ans(4)^2*As(1)*As(2) - 1536*Ans(5)*Ans(1)*Ans(2)*Ans(4)*As(1)^2 + 1024*Ans(1)*Ans(3)^3*As(1)*As(3) - 512*Ans(1)*Ans(3)^2*Ans(4)*As(1)*As(2) + 1024*Ans(5)*Ans(1)*Ans(3)^2*As(1)^2 - 432*Ans(2)^4*As(3)^2 + 288*Ans(2)^3*Ans(3)*As(2)*As(3) + 288*Ans(2)^3*Ans(4)*As(1)*As(3) - 144*Ans(2)^3*Ans(4)*As(2)^2 - 192*Ans(2)^2*Ans(3)^2*As(1)*As(3) + 96*Ans(2)^2*Ans(3)*Ans(4)*As(1)*As(2) - 48*Ans(2)^2*Ans(4)^2*As(1)^2;
k3 = - 4096*Ans(1)^3*As(3)^3 + 3072*Ans(1)^2*Ans(2)*As(2)*As(3)^2 + 4096*Ans(1)^2*Ans(3)*As(1)*As(3)^2 - 2048*Ans(1)^2*Ans(3)*As(2)^2*As(3) - 1024*Ans(5)*Ans(1)^2*As(1)*As(2)^2 - 3072*Ans(4)*Ans(1)^2*As(1)*As(2)*As(3) + 1280*Ans(4)*Ans(1)^2*As(2)^3 - 2304*Ans(1)*Ans(2)^2*As(1)*As(3)^2 + 1536*Ans(1)*Ans(2)*Ans(3)*As(1)*As(2)*As(3) - 128*Ans(1)*Ans(2)*Ans(3)*As(2)^3 + 1024*Ans(5)*Ans(1)*Ans(2)*As(1)^2*As(2) + 1536*Ans(4)*Ans(1)*Ans(2)*As(1)^2*As(3) - 896*Ans(4)*Ans(1)*Ans(2)*As(1)*As(2)^2 - 1024*Ans(1)*Ans(3)^2*As(1)^2*As(3) + 256*Ans(1)*Ans(3)^2*As(1)*As(2)^2 + 32*Ans(2)^3*As(2)^3 - 64*Ans(2)^2*Ans(3)*As(1)*As(2)^2 - 256*Ans(5)*Ans(2)^2*As(1)^3 + 128*Ans(4)*Ans(2)^2*As(1)^2*As(2);
k1 = 1024*As(3)*Ans(1)^2*As(1)*As(2)^2 - 256*Ans(1)^2*As(2)^4 - 1024*As(3)*Ans(1)*Ans(2)*As(1)^2*As(2) + 256*Ans(1)*Ans(2)*As(1)*As(2)^3 + 256*As(3)*Ans(2)^2*As(1)^3 - 64*Ans(2)^2*As(1)^2*As(2)^2;



end

