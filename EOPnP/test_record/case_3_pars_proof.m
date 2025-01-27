syms a1 a2 a3 k1 k2 real
syms v11 v12 v13 v14 v15 v16 v17 v18 v19 real;
syms v21 v22 v23 v24 v25 v26 v27 v28 v29 real;
syms v31 v32 v33 v34 v35 v36 v37 v38 v39 real;
v1 = [v11;v12;v13;v14;v15;v16;v17;v18;v19];
v2 = [v21;v22;v23;v24;v25;v26;v27;v28;v29];
v3 = [v31;v32;v33;v34;v35;v36;v37;v38;v39];

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

%求矩阵B，并证明用dot运算可以组成矩阵B的元素
r1 = a1*v1_3 + a2*v2_3 + a3*v3_3;
r2 = a1*v1_6 + a2*v2_6 + a3*v3_6;
r3 = a1*v1_9 + a2*v2_9 + a3*v3_9;

g1 = dot(r1,r1) - dot(r2,r2);
g1 = collect(g1,[a1 a2 a3]);
% (v11^2 + v12^2 + v13^2 - v14^2 - v15^2 - v16^2)*a1^2 + (2*v11*v21 + 2*v12*v22 + 2*v13*v23 - 2*v14*v24 - 2*v15*v25 - 2*v16*v26)*a1*a2
% + (2*v11*v31 + 2*v12*v32 + 2*v13*v33 - 2*v14*v34 - 2*v15*v35 - 2*v16*v36)*a1*a3 + (v21^2 + v22^2 + v23^2 - v24^2 - v25^2 - v26^2)*a2^2
% + (2*v21*v31 + 2*v22*v32 + 2*v23*v33 - 2*v24*v34 - 2*v25*v35 - 2*v26*v36)*a2*a3 + (v31^2 + v32^2 + v33^2 - v34^2 - v35^2 - v36^2)*a3^2

if( dot(v1_3,v1_3) - dot(v1_6,v1_6) == v11^2 + v12^2 + v13^2 - v14^2 - v15^2 - v16^2) 
    disp('1equal');  
end
% ( dot(v1_3,v1_3) - dot(v1_6,v1_6) ) * a1^2

if( 2*( dot(v1_3,v2_3) - dot(v1_6,v2_6) ) == 2*v11*v21 + 2*v12*v22 + 2*v13*v23 - 2*v14*v24 - 2*v15*v25 - 2*v16*v26) 
    disp('2equal');  
end
% 2*( dot(v1_3,v2_3) - dot(v1_6,v2_6) ) * a1 * a2

if( 2*( dot(v1_3,v3_3) - dot(v1_6,v3_6) ) == 2*v11*v31 + 2*v12*v32 + 2*v13*v33 - 2*v14*v34 - 2*v15*v35 - 2*v16*v36) 
    disp('3equal');  
end
% 2*( dot(v1_3,v3_3) - dot(v1_6,v3_6) ) * a1 * a3

if( dot(v2_3,v2_3) - dot(v2_6,v2_6) == v21^2 + v22^2 + v23^2 - v24^2 - v25^2 - v26^2) 
    disp('4equal');  
end
% ( dot(v2_3,v2_3) - dot(v2_6,v2_6) ) * a2^2

if( 2*( dot(v2_3,v3_3) - dot(v2_6,v3_6) ) == 2*v21*v31 + 2*v22*v32 + 2*v23*v33 - 2*v24*v34 - 2*v25*v35 - 2*v26*v36) 
    disp('5equal');  
end
% 2*( dot(v2_3,v3_3) - dot(v2_6,v3_6) ) * a2 * a3

if( dot(v3_3,v3_3) - dot(v3_6,v3_6) == v31^2 + v32^2 + v33^2 - v34^2 - v35^2 - v36^2) 
    disp('6equal');  
end
% ( dot(v3_3,v3_3) - dot(v3_6,v3_6) ) * a3^2

%证得 g1 = dot(r1,r1) - dot(r2,r2)
%        = ( dot(v1_3,v1_3) - dot(v1_6,v1_6) ) * a1^2 + 2*( dot(v1_3,v2_3)
%        - dot(v1_6,v2_6) ) * a1 * a2 + 2*( dot(v1_3,v3_3) - dot(v1_6,v3_6)
%        ) * a1 * a3 + ( dot(v2_3,v2_3) - dot(v2_6,v2_6) ) * a2^2 + 2*(
%        dot(v2_3,v3_3) - dot(v2_6,v3_6) ) * a2 * a3 + ( dot(v3_3,v3_3) -
%        dot(v3_6,v3_6) ) * a3^2 

%同理 g2 = dot(r1,r1) - dot(r3,r3)
%        = ( dot(v1_3,v1_3) - dot(v1_9,v1_9) ) * a1^2 + 2*( dot(v1_3,v2_3)
%        - dot(v1_9,v2_9) ) * a1 * a2 + 2*( dot(v1_3,v3_3) - dot(v1_9,v3_9)
%        ) * a1 * a3 + ( dot(v2_3,v2_3) - dot(v2_9,v2_9) ) * a2^2 + 2*(
%        dot(v2_3,v3_3) - dot(v2_9,v3_9) ) * a2 * a3 + ( dot(v3_3,v3_3) -
%        dot(v3_9,v3_9) ) * a3^2 

g2_true = dot(r1,r1) - dot(r3,r3);
g2_e = ( dot(v1_3,v1_3) - dot(v1_9,v1_9) ) * a1^2 + 2*( dot(v1_3,v2_3) - dot(v1_9,v2_9) ) * a1 * a2 + 2*( dot(v1_3,v3_3) - dot(v1_9,v3_9) ) * a1 * a3 + ( dot(v2_3,v2_3) - dot(v2_9,v2_9) ) * a2^2 + 2*( dot(v2_3,v3_3) - dot(v2_9,v3_9) ) * a2 * a3 + ( dot(v3_3,v3_3) - dot(v3_9,v3_9) ) * a3^2; 
if( simplify(g2_true - g2_e) == 0 )
    disp('7equal');
end
%证得f2表达式

f4 = dot(r1,r2);
f4 = collect(f4,[a1 a2 a3]);
%(v11*v14 + v12*v15 + v13*v16)*a1^2 + (v11*v24 + v14*v21 + v12*v25 + v15*v22 + v13*v26 + v16*v23)*a1*a2
% + (v11*v34 + v14*v31 + v12*v35 + v15*v32 + v13*v36 + v16*v33)*a1*a3 + (v21*v24 + v22*v25 + v23*v26)*a2^2
% + (v21*v34 + v24*v31 + v22*v35 + v25*v32 + v23*v36 + v26*v33)*a2*a3 + (v31*v34 + v32*v35 + v33*v36)*a3^2

f4_e = dot(v1_3,v1_6)*a1^2 + ( dot(v1_3,v2_6) + dot(v1_6,v2_3) )*a1*a2 + ( dot(v1_3,v3_6) + dot(v1_6,v3_3) )*a1*a3 + ( dot(v2_3,v3_6) + dot(v2_6,v3_3) )*a2*a3 + dot(v2_3,v2_6)*a2^2 + dot(v3_3,v3_6)*a3^2;
if( simplify(f4 - f4_e) == 0 )
    disp('8equal');
end
%证得 f4 = dot(v1_3,v1_6)*a1^2 + ( dot(v1_3,v2_6) + dot(v1_6,v2_3) )*a1*a2
%           + ( dot(v1_3,v3_6) + dot(v1_6,v3_3) )*a1*a3 + ( dot(v2_3,v3_6)
%           + dot(v2_6,v3_3) )*a2*a3 + dot(v2_3,v2_6)*a2^2 + dot(v3_3,v3_6)*a3^2;

f5 = dot(r1,r3);
f5 = collect(f5,[a1 a2 a3]);
f5_e = subs(f4,[v1_6 v2_6 v3_6],[v1_9 v2_9 v3_9]);
if( simplify(f5 - f5_e) == 0 )
    disp('9equal');
end
%证得 f5 = dot(v1_3,v1_9)*a1^2 + ( dot(v1_3,v2_9) + dot(v1_9,v2_3) )*a1*a2
%           + ( dot(v1_3,v3_9) + dot(v1_9,v3_3) )*a1*a3 + ( dot(v2_3,v3_9)
%           + dot(v2_9,v3_3) )*a2*a3 + dot(v2_3,v2_9)*a2^2 + dot(v3_3,v3_9)*a3^2;

f6 = dot(r2,r3);
f6 = collect(f6,[a1 a2 a3]);
f6_e = subs(f5,[v1_3 v2_3 v3_3],[v1_6 v2_6 v3_6]);
if( simplify(f6 - f6_e) == 0 )
    disp('10equal');
end
%证得 f6 = dot(v1_6,v1_9)*a1^2 + ( dot(v1_6,v2_9) + dot(v1_9,v2_6) )*a1*a2
%           + ( dot(v1_6,v3_9) + dot(v1_9,v3_6) )*a1*a3 + ( dot(v2_6,v3_9)
%           + dot(v2_9,v3_6) )*a2*a3 + dot(v2_6,v2_9)*a2^2 + dot(v3_6,v3_9)*a3^2;

%计算B矩阵
alpha = [ a1^2; a2^2; a3^2; a1*a2; a1*a3; a2*a3 ];
B1 = [   dot(v1_3,v1_3) - dot(v1_6,v1_6)     dot(v2_3,v2_3) - dot(v2_6,v2_6)     dot(v3_3,v3_3) - dot(v3_6,v3_6)     2*( dot(v1_3,v2_3) - dot(v1_6,v2_6) )    2*( dot(v1_3,v3_3) - dot(v1_6,v3_6) )    2*( dot(v2_3,v3_3) - dot(v2_6,v3_6) )  ];
if( simplify(g1 - B1*alpha) == 0 )
    disp('11equal');
end       

B2 = [   dot(v1_3,v1_3) - dot(v1_9,v1_9)     dot(v2_3,v2_3) - dot(v2_9,v2_9)     dot(v3_3,v3_3) - dot(v3_9,v3_9)     2*( dot(v1_3,v2_3) - dot(v1_9,v2_9) )    2*( dot(v1_3,v3_3) - dot(v1_9,v3_9) )    2*( dot(v2_3,v3_3) - dot(v2_9,v3_9) )  ];
g2 = dot(r1,r1) - dot(r3,r3);
if( simplify(g2 - B2*alpha) == 0 )
    disp('12equal');
end    

B3 = [	dot(v1_3,v1_6)  	dot(v2_3,v2_6)  	dot(v3_3,v3_6)    dot(v1_3,v2_6) + dot(v1_6,v2_3) 	 dot(v1_3,v3_6) + dot(v1_6,v3_3) 	 dot(v2_3,v3_6) + dot(v2_6,v3_3) ];		
if( simplify(f4 - B3*alpha) == 0 )
    disp('13equal');
end 

B4 = [	dot(v1_3,v1_9)  	dot(v2_3,v2_9)  	dot(v3_3,v3_9)    dot(v1_3,v2_9) + dot(v1_9,v2_3) 	 dot(v1_3,v3_9) + dot(v1_9,v3_3) 	 dot(v2_3,v3_9) + dot(v2_9,v3_3) ];		
f5 = dot(r1,r3);
if( simplify(f5 - B4*alpha) == 0 )
    disp('14equal');
end 

B5 = [	dot(v1_6,v1_9)  	dot(v2_6,v2_9)  	dot(v3_6,v3_9)    dot(v1_6,v2_9) + dot(v1_9,v2_6) 	 dot(v1_6,v3_9) + dot(v1_9,v3_6) 	 dot(v2_6,v3_9) + dot(v2_9,v3_6) ];		
f6 = dot(r2,r3);
if( simplify(f6 - B5*alpha) == 0 )
    disp('15equal');
end 

B = [      dot(v1_3,v1_3) - dot(v1_6,v1_6)     dot(v2_3,v2_3) - dot(v2_6,v2_6)     dot(v3_3,v3_3) - dot(v3_6,v3_6)     2*( dot(v1_3,v2_3) - dot(v1_6,v2_6) )    2*( dot(v1_3,v3_3) - dot(v1_6,v3_6) )      2*( dot(v2_3,v3_3) - dot(v2_6,v3_6) )   ;
	       dot(v1_3,v1_3) - dot(v1_9,v1_9)     dot(v2_3,v2_3) - dot(v2_9,v2_9)     dot(v3_3,v3_3) - dot(v3_9,v3_9)     2*( dot(v1_3,v2_3) - dot(v1_9,v2_9) )    2*( dot(v1_3,v3_3) - dot(v1_9,v3_9) )      2*( dot(v2_3,v3_3) - dot(v2_9,v3_9) )   ;
	       dot(v1_3,v1_6)                      dot(v2_3,v2_6)                      dot(v3_3,v3_6)                      dot(v1_3,v2_6) + dot(v1_6,v2_3)         dot(v1_3,v3_6) + dot(v1_6,v3_3)             dot(v2_3,v3_6) + dot(v2_6,v3_3)         ;		
	       dot(v1_3,v1_9)                      dot(v2_3,v2_9)                      dot(v3_3,v3_9)                      dot(v1_3,v2_9) + dot(v1_9,v2_3)         dot(v1_3,v3_9) + dot(v1_9,v3_3)             dot(v2_3,v3_9) + dot(v2_9,v3_3)         ;		
	       dot(v1_6,v1_9)                      dot(v2_6,v2_9)                      dot(v3_6,v3_9)                      dot(v1_6,v2_9) + dot(v1_9,v2_6)         dot(v1_6,v3_9) + dot(v1_9,v3_6)             dot(v2_6,v3_9) + dot(v2_9,v3_6)         ];		



%将alpha向量换成k向量
g1_k = subs(g1,[a2 a3],[k1*a1 k2*a1]);
g1_k = g1_k/(a1^2);
k = subs(alpha,[a2 a3],[k1*a1 k2*a1]);
k = k./(a1^2);
% k =[ 1; k1^2; k2^2; k1; k2; k1*k2 ]
if( simplify(g1_k - B1*k) == 0 )
    disp('16equal');
end 
%证得 k =[ 1; k1^2; k2^2; k1; k2; k1*k2 ] 
%注意 论文中k元素顺序出错（若B保持不变）

%求C_k1k2表达式
%求(B1*k)^2形式
syms b1 b2 b3 b4 b5 b6
b = [ b1 b2 b3 b4 b5 b6 ];
bk2 = (b*k)^2;
bk2 = collect(bk2,[k1 k2]);
[b_new k_new] = coeffs(bk2,[k1 k2]);
k_new = k_new';
% k_new = [ k1^4; k1^3*k2; k1^3; k1^2*k2^2; k1^2*k2; k1^2; k1*k2^3; k1*k2^2; k1*k2; k1; k2^4; k2^3; k2^2; k2; 1]


%根据上述形式求B_new和Bns
% B1_new = subs(b_new,b,B1);
% B2_new = subs(b_new,b,B2);
% B3_new = subs(b_new,b,B3);
% B4_new = subs(b_new,b,B4);
% B5_new = subs(b_new,b,B5);
% B_new = [B1_new;B2_new;B3_new;B4_new;B5_new];   %5*15矩阵
sym('B_new',[5,15],'real');
for i = 1 : 5
    B_new(i,1) = B(i,2)^2;
    B_new(i,2) = 2*B(i,2)*B(i,6);
    B_new(i,3) = 2*B(i,2)*B(i,4);
    B_new(i,4) = B(i,6)^2 + 2*B(i,2)*B(i,3); 
    B_new(i,5) = 2*B(i,2)*B(i,5) + 2*B(i,4)*B(i,6);
    B_new(i,6) = B(i,4)^2 + 2*B(i,1)*B(i,2);
    B_new(i,7) = 2*B(i,3)*B(i,6);
    B_new(i,8) = 2*B(i,3)*B(i,4) + 2*B(i,5)*B(i,6);
    B_new(i,9) = 2*B(i,1)*B(i,6) + 2*B(i,4)*B(i,5);
    B_new(i,10) = 2*B(i,1)*B(i,4);
    B_new(i,11) = B(i,3)^2;
    B_new(i,12) = 2*B(i,3)*B(i,5);
    B_new(i,13) = B(i,5)^2 + 2*B(i,1)*B(i,3);
    B_new(i,14) = 2*B(i,1)*B(i,5);
    B_new(i,15) = B(i,1)^2;
end
Bns = B_new(1,:) + B_new(2,:) + B_new(3,:) + B_new(4,:) + B_new(5,:);   %1*15矩阵

%证明
C_k1k2 = (B1*k)^2 + (B2*k)^2 + (B3*k)^2 + (B4*k)^2 + (B5*k)^2;
t1 = diff(C_k1k2,k1);
t1_e = Bns * diff(k_new,k1);
if( simplify(t1 - t1_e) == 0 )
    disp('17equal');
end 
%证得 B_new、Bns正确 且 t1 = Bns * diff(k_new,k1);
t2 = diff(C_k1k2,k2);
t2_e = Bns * diff(k_new,k2);
if( simplify(t2 - t2_e) == 0 )
    disp('18equal');
end 
%证得 B_new、Bns正确 且 t2 = Bns * diff(k_new,k2);        

%证明c11等系数形式
%先求真值
c1 = coeffs(t1,k1);
c2 = coeffs(t2,k1);

c11 = c1(1,4);
c12 = c1(1,3);
c13 = c1(1,2);
c14 = c1(1,1);

c21 = c2(1,4);
c22 = c2(1,3);
c23 = c2(1,2);
c24 = c2(1,1);

%Bns形式表达c11
coe11 = 4*Bns(1);
coe12 = 3*k2*Bns(2) + 3*Bns(3);
coe13 = 2*k2^2*Bns(4) + 2*k2*Bns(5) + 2*Bns(6);
coe14 = k2^3*Bns(7) + k2^2*Bns(8) + k2*Bns(9) + Bns(10);

f7 = [ coe11 coe12 coe13 coe14 ] - [ c11 c12 c13 c14 ];
if( simplify(f7) == 0 )
    disp('19equal');
end 
%证得 

coe21 = Bns(2);
coe22 = 2*k2*Bns(4) + Bns(5);
coe23 = 3*k2^2*Bns(7) + 2*k2*Bns(8) + Bns(9);
coe24 = 4*k2^3*Bns(11) + 3*k2^2*Bns(12) + 2*k2*Bns(13) + Bns(14);
  
f8 = [ coe21 coe22 coe23 coe24 ] - [ c21 c22 c23 c24 ];
if( simplify(f8) == 0 )
    disp('20equal');
end 
%证得 
 


%上面所证得的Bns为已知部分，可以由输入数据运算得出
%将Bns作为符号变量固定，求最终系数k9、k7等
syms Bnss1 Bnss2 Bnss3 Bnss4 Bnss5 Bnss6 Bnss7 Bnss8 Bnss9 Bnss10 Bnss11 Bnss12 Bnss13 Bnss14 Bnss15 real
Bnss = [ Bnss1 Bnss2 Bnss3 Bnss4 Bnss5 Bnss6 Bnss7 Bnss8 Bnss9 Bnss10 Bnss11 Bnss12 Bnss13 Bnss14 Bnss15 ];

coe11 = 4*Bnss(1);
coe12 = 3*k2*Bnss(2) + 3*Bnss(3);
coe13 = 2*k2^2*Bnss(4) + 2*k2*Bnss(5) + 2*Bnss(6);
coe14 = k2^3*Bnss(7) + k2^2*Bnss(8) + k2*Bnss(9) + Bnss(10);

coe21 = Bnss(2);
coe22 = 2*k2*Bnss(4) + Bnss(5);
coe23 = 3*k2^2*Bnss(7) + 2*k2*Bnss(8) + Bnss(9);
coe24 = 4*k2^3*Bnss(11) + 3*k2^2*Bnss(12) + 2*k2*Bnss(13) + Bnss(14);

mat1 = [coe11    0          0           coe21   0       0     ;
        coe12    coe11      0           coe22   coe21   0     ;
        coe13    coe12      coe11       coe23   coe22   coe21 ;
        coe14    coe13      coe12       coe24   coe23   coe22 ;
          0       coe14     coe13       0       coe24   coe23 ;
          0       0         coe14       0       0       coe24  ];

f_k = det(mat1);
[c,s] = coeffs(f_k,k2);
k9 = c(1);
k8 = c(2);
k7 = c(3);
k6 = c(4);
k5 = c(5);
k4 = c(6);
k3 = c(7);
kk2 = c(8);
kk1 = c(9);
k0 = c(10);

%求代回f1形式
f1 = dot(r1,r1);
f1 = subs(f1,[a2 a3],[k1*a1 k2*a1]);
f1 = collect(f1,a1);
[ca1,sa1] = coeffs(f1,a1)

