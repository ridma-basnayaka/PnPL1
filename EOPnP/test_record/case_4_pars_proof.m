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

% r1 = a1*v1_3 + a2*v2_3 + a3*v3_3 + a4*v4_3;
% r2 = a1*v1_6 + a2*v2_6 + a3*v3_6 + a4*v4_6;
% r3 = a1*v1_9 + a2*v2_9 + a3*v3_9 + a4*v4_9;

r1 = v1_3 + k1*v2_3 + k2*v3_3 + k3*v4_3;
r2 = v1_6 + k1*v2_6 + k2*v3_6 + k3*v4_6;
r3 = v1_9 + k1*v2_9 + k2*v3_9 + k3*v4_9;

f4k = dot(r1,r2);
[cf4k,sf4k] = coeffs(f4k,[k1 k2 k3]);
f5k = dot(r1,r3);
f6k = dot(r2,r3);

%验证C矩阵和k矩阵形式
C = [   dot(v2_3,v2_6)  dot(v3_3,v3_6)  dot(v4_3,v4_6)  dot(v2_3,v3_6)+dot(v2_6,v3_3)   dot(v2_3,v4_6)+dot(v2_6,v4_3)   dot(v3_3,v4_6)+dot(v3_6,v4_3)   dot(v1_3,v2_6)+dot(v1_6,v2_3)   dot(v1_3,v3_6)+dot(v1_6,v3_3)   dot(v1_3,v4_6)+dot(v1_6,v4_3)   dot(v1_3,v1_6);
        dot(v2_3,v2_9)  dot(v3_3,v3_9)  dot(v4_3,v4_9)  dot(v2_3,v3_9)+dot(v2_9,v3_3)   dot(v2_3,v4_9)+dot(v2_9,v4_3)   dot(v3_3,v4_9)+dot(v3_9,v4_3)   dot(v1_3,v2_9)+dot(v1_9,v2_3)   dot(v1_3,v3_9)+dot(v1_9,v3_3)   dot(v1_3,v4_9)+dot(v1_9,v4_3)   dot(v1_3,v1_9);
        dot(v2_6,v2_9)  dot(v3_6,v3_9)  dot(v4_6,v4_9)  dot(v2_6,v3_9)+dot(v2_9,v3_6)   dot(v2_6,v4_9)+dot(v2_9,v4_6)   dot(v3_6,v4_9)+dot(v3_9,v4_6)   dot(v1_6,v2_9)+dot(v1_9,v2_6)   dot(v1_6,v3_9)+dot(v1_9,v3_6)   dot(v1_6,v4_9)+dot(v1_9,v4_6)   dot(v1_6,v1_9)];

k = [ k1^2; k2^2; k3^2; k1*k2; k1*k3; k2*k3; k1; k2; k3; 1];
        
if( simplify( C*k - [ f4k; f5k; f6k ] ) == 0)
    disp('equal');
end

%固定C，求M矩阵和N_k矩阵形式
syms C11 C12 C13 C14 C15 C16 C17 C18 C19 C110 real
syms C21 C22 C23 C24 C25 C26 C27 C28 C29 C210 real
syms C31 C32 C33 C34 C35 C36 C37 C38 C39 C310 real

Cv = [  C11 C12 C13 C14 C15 C16 C17 C18 C19 C110;
        C21 C22 C23 C24 C25 C26 C27 C28 C29 C210;
        C31 C32 C33 C34 C35 C36 C37 C38 C39 C310 ];

f = Cv*k;

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
 
Nk = (-1) * [f4_k3k2(1,2) f5_k3k2(1,2) f6_k3k2(1,2);
              f4_k2k3(1,2) f5_k2k3(1,2) f6_k2k3(1,2);
              f4_k3k2(1,1) f5_k3k2(1,1) f6_k3k2(1,1)]';
          
% M =
%  
% [ C12, C13, C16]
% [ C22, C23, C26]
% [ C32, C33, C36]

% Nk =
%  
% [ - C18 - C14*k1, - C19 - C15*k1, - C11*k1^2 - C17*k1 - C110]
% [ - C28 - C24*k1, - C29 - C25*k1, - C21*k1^2 - C27*k1 - C210]
% [ - C38 - C34*k1, - C39 - C35*k1, - C31*k1^2 - C37*k1 - C310]
          

%M、pinv(M)已知，固定pinv(M)
%Nk可以分为三部分相加（笔记），和Mp相乘后，得到的Tk也可以由三部分相加表示，将其中已知部分用符号变量固定
syms X11 X12 X13 X21 X22 X23 X31 X32 X33 real
X = [X11 X12 X13; X21 X22 X23; X31 X32 X33];
syms Y11 Y12 Y13 Y21 Y22 Y23 Y31 Y32 Y33 real
Y = [Y11 Y12 Y13; Y21 Y22 Y23; Y31 Y32 Y33];
syms Z1 Z2 Z3 real
Z = [Z1;Z2;Z3];

Tk = X + k1.*Y + Z*[0 0 k1^2];

K = Tk * [k2;k3;1];
k2_2 = K(1,:);
k3_2 = K(2,:);
k2_k3 = K(3,:);

M_k(1,1) = k2_k3*Tk(1,1) + k3_2*Tk(1,2) + k3*Tk(1,3) - k2_2*Tk(3,1) - k2_k3*Tk(3,2) -k2*Tk(3,3);
M_k(2,1) = k2_k3*Tk(3,1) + k3_2*Tk(3,2) + k3*Tk(3,3) - k2_2*Tk(2,1) - k2_k3*Tk(2,2) -k2*Tk(2,3);
M_k(3,1) = k2_2*Tk(3,1)^2 + k2_k3*( 2*Tk(3,1)*Tk(3,2) ) + k3_2*Tk(3,2)^2 + 2*k2*Tk(3,1)*Tk(3,3) + 2*k3*Tk(3,2)*Tk(3,3) + Tk(3,3)^2 - ( k2_2*Tk(1,1)*Tk(2,1) + k2_k3*( Tk(1,1)*Tk(1,2) + Tk(1,2)*Tk(2,1) ) + k3_2*Tk(1,2)*Tk(2,2) + k2*( Tk(1,1)*Tk(2,3) + Tk(2,1)*Tk(2,3) )  + k3*( Tk(1,2)*Tk(2,3) + Tk(2,2)*Tk(1,3) ) + Tk(1,3)*Tk(2,3) );

%提取系数矩阵
[s1,~] = coeffs(M_k(1),[k2 k3]);
[s2,~] = coeffs(M_k(2),[k2 k3]);
[s3,~] = coeffs(M_k(3),[k2 k3]);
M_k = [s1;s2;s3];
detMk = M_k(1,1)*M_k(2,2)*M_k(3,3) + M_k(1,2)*M_k(2,3)*M_k(3,1) + M_k(1,3)*M_k(2,1)*M_k(3,2) - M_k(1,3)*M_k(2,2)*M_k(3,1) - M_k(1,1)*M_k(2,3)*M_k(3,2) - M_k(3,3)*M_k(2,1)*M_k(1,2);
[cdm,sdm] = coeffs(detMk,k1);

%求a1代价方程系数
r1 = a1*v1_3 + k1*a1*v2_3 + k2*a1*v3_3 + k3*a1*v4_3;
r2 = a1*v1_6 + k1*a1*v2_6 + k2*a1*v3_6 + k3*a1*v4_6;
r3 = a1*v1_9 + k1*a1*v2_9 + k2*a1*v3_9 + k3*a1*v4_9;
fa1 = (dot(r1,r1) - 1)^2 + (dot(r2,r2) - 1)^2 + (dot(r3,r3) - 1)^2;
fa1 = diff(fa1,a1);
[cfa,sfa] = coeffs(fa1,a1)















