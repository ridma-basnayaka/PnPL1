syms a1 a2 real
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

%求矩阵A，并证明用dot运算可以组成矩阵A的元素
r1 = a1*v1_3 + a2*v2_3;
r2 = a1*v1_6 + a2*v2_6;
r3 = a1*v1_9 + a2*v2_9;

dr1r1 = dot(r1,r1);
dr1r1 = collect(dr1r1,[a1 a2]); 
% result: (v11^2 + v12^2 + v13^2)*a1^2 + (2*v11*v21 + 2*v12*v22 + 2*v13*v23)*a1*a2 + (v21^2 + v22^2 + v23^2)*a2^2
if(dot(v1_3,v1_3) == v11^2 + v12^2 + v13^2) 
    disp('1equal');  
end
if(2*dot(v1_3,v2_3) == 2*v11*v21 + 2*v12*v22 + 2*v13*v23) 
    disp('2equal');  
end
if(dot(v2_3,v2_3) == v21^2 + v22^2 + v23^2) 
    disp('3equal');  
end
% 证明 dot(r1,r1) = dot(v1_3,v1_3)*a1^2 + 2*dot(v1_3,v2_3)*a1*a2 + dot(v2_3,v2_3)*a2^2
% 同理 dot(r2,r2) = dot(v1_6,v1_6)*a1^2 + 2*dot(v1_6,v2_6)*a1*a2 + dot(v2_6,v2_6)*a2^2
% 同理 dot(r3,r3) = dot(v1_9,v1_9)*a1^2 + 2*dot(v1_9,v2_9)*a1*a2 + dot(v2_9,v2_9)*a2^2

dr1r2 = dot(r1,r2);
dr1r2 = collect(dr1r2,[a1 a2]);
% result: (v11*v14 + v12*v15 + v13*v16)*a1^2 + (v11*v24 + v14*v21 + v12*v25 + v15*v22 + v13*v26 + v16*v23)*a1*a2 + (v21*v24 + v22*v25 + v23*v26)*a2^2
if(dot(v1_3,v1_6) == v11*v14 + v12*v15 + v13*v16) 
    disp('4equal');  
end
if(dot(v1_3,v2_6) + dot(v2_3,v1_6) == v11*v24 + v14*v21 + v12*v25 + v15*v22 + v13*v26 + v16*v23) 
    disp('5equal');  
end
if(dot(v2_6,v2_3) == v21*v24 + v22*v25 + v23*v26) 
    disp('6equal');  
end
% 证明 dot(r1,r2) = dot(v1_3,v1_6)*a1^2 + ( dot(v1_3,v2_6) + dot(v2_3,v1_6) )*a1*a2 + dot(v2_3,v2_6)*a2^2
% 同理 dot(r1,r3) = dot(v1_3,v1_9)*a1^2 + ( dot(v1_3,v2_9) + dot(v2_3,v1_9) )*a1*a2 + dot(v2_3,v2_9)*a2^2
% 同理 dot(r2,r3) = dot(v1_6,v1_9)*a1^2 + ( dot(v1_6,v2_9) + dot(v2_6,v1_9) )*a1*a2 + dot(v2_6,v2_9)*a2^2

% 证得 A = 
%       [ dot(v1_3,v1_3)    2 * dot(v1_3,v2_3)                      dot(v2_3,v2_3);	
%         dot(v1_6,v1_6)    2 * dot(v1_6,v2_6)                      dot(v2_6,v2_6);	
%         dot(v1_9,v1_9)    2 * dot(v1_9,v2_9)                      dot(v2_9,v2_9);	
%         dot(v1_3,v1_6)	( dot(v1_3,v2_6) + dot(v2_3,v1_6) )     dot(v2_3,v2_6);
%         dot(v1_3,v1_9)	( dot(v1_3,v2_9) + dot(v2_3,v1_9) )     dot(v2_3,v2_9);
%         dot(v1_6,v1_9)	( dot(v1_6,v2_9) + dot(v2_6,v1_9) )     dot(v2_6,v2_9) ]
%
%   alpha = [a1^2;a1*a2;a2^2]
%   b = [1;1;1;0;0;0]
%   A * alpha = b


f1 = ( dot(r1,r1) - 1 )^2;
f2 = dot(r1,r1)^2 - 2*dot(r1,r1) + 1;
f1 = collect(f1,[a1 a2]);
f2 = collect(f2,[a1 a2]);
if( f1 - f2 == 0 ) 
    disp('7equal');  
end
%证得  ( dot(r1,r1) - 1 )^2 = dot(r1,r1)^2 - 2*dot(r1,r1) + 1;
% ( dot(r2,r2) - 1 )^2 同理
% ( dot(r3,r3) - 1 )^2 同理



A = [   dot(v1_3,v1_3)    2 * dot(v1_3,v2_3)                        dot(v2_3,v2_3);	
        dot(v1_6,v1_6)    2 * dot(v1_6,v2_6)                        dot(v2_6,v2_6);	
        dot(v1_9,v1_9)    2 * dot(v1_9,v2_9)                        dot(v2_9,v2_9);	
        dot(v1_3,v1_6)	( dot(v1_3,v2_6) + dot(v2_3,v1_6) )         dot(v2_3,v2_6);
        dot(v1_3,v1_9)	( dot(v1_3,v2_9) + dot(v2_3,v1_9) )         dot(v2_3,v2_9);
        dot(v1_6,v1_9)	( dot(v1_6,v2_9) + dot(v2_6,v1_9) )         dot(v2_6,v2_9) ];
alpha = [a1^2;a1*a2;a2^2];
f3 = ( A(1,:)*alpha )^2 - 2*( A(1,:)*alpha ) + 1;
f3 = collect(f3,[a1 a2]);
if( f1 - f3 == 0 ) 
    disp('8equal');  
end
%证得  ( dot(r1,r1) - 1 )^2 = ( A(1,:)*alpha )^2 - 2*( A(1,:)*alpha ) + 1;
%同理  ( dot(r2,r2) - 1 )^2 = ( A(2,:)*alpha )^2 - 2*( A(2,:)*alpha ) + 1;
%      ( dot(r3,r3) - 1 )^2 = ( A(3,:)*alpha )^2 - 2*( A(3,:)*alpha ) + 1;

syms A11 A12 A13
A1 = [A11 A12 A13];
f4 = ( A1*alpha )^2;
f4 = collect(f4,[a1 a2]);
% result: A11^2*a1^4 + 2*A11*A12*a1^3*a2 + (A12^2 + 2*A11*A13)*a1^2*a2^2 + 2*A12*A13*a1*a2^3 + A13^2*a2^4
% 所以 alpha_new = [a1^4 ; a1^3*a2 ; a1^2*a2^2 ; a1*a2^3 ; a2^4]
% A1_new = [ A11^2   2*A11*A12   A12^2 + 2*A11*A13   2*A12*A13   A13^2]

C_a1a2 = ( dot(r1,r1) - 1 )^2 + ( dot(r2,r2) - 1 )^2 + ( dot(r3,r3) - 1 )^2 + dot(r1,r2)^2 + dot(r1,r3)^2 + dot(r2,r3)^2;
C_a1a2 = collect(C_a1a2,[a1 a2]);
f5 = ( A(1,:)*alpha )^2 + ( A(2,:)*alpha )^2 + ( A(3,:)*alpha )^2 + ( A(4,:)*alpha )^2 + ( A(5,:)*alpha )^2 + ( A(6,:)*alpha )^2 + 3 - 2*( A(1,:)*alpha + A(2,:)*alpha + A(3,:)*alpha );
f5 = collect(f5,[a1 a2]);
f6 = alpha'*A'*A*alpha + 3 - 2*( A(1,:)*alpha + A(2,:)*alpha + A(3,:)*alpha ); 
f6 = collect(f6,[a1 a2]);
if( C_a1a2 - f5 == 0 ) 
    disp('9equal');  
end
if( C_a1a2 - f6 == 0 ) 
    disp('10equal');  
end
%证得 C_a1a2 = f6 = f5

%根据上面对( A1*alpha )^2的转化
sym('A_new',[6,5],'real');
for i = 1 : 6
    A_new(i,1) = A(i,1)^2;
    A_new(i,2) = 2*A(i,1)*A(i,2);
    A_new(i,3) = A(i,2)^2 + 2*A(i,1)*A(i,3);
    A_new(i,4) = 2*A(i,2)*A(i,3); 
    A_new(i,5) = A(i,3)^2;
end
alpha_new = [a1^4 ; a1^3*a2 ; a1^2*a2^2 ; a1*a2^3 ; a2^4];
f7 = ( A_new(1,:) + A_new(2,:) + A_new(3,:) + A_new(4,:) + A_new(5,:) + A_new(6,:) )* alpha_new; 
f7 = collect(f7,[a1 a2]);
f8 = ( A(1,:)*alpha )^2 + ( A(2,:)*alpha )^2 + ( A(3,:)*alpha )^2 + ( A(4,:)*alpha )^2 + ( A(5,:)*alpha )^2 + ( A(6,:)*alpha )^2;
f8 = collect(f8,[a1 a2]);
if ( simplify(f7 - f8) == 0 )
    disp('11equal');  
end
%证得 f7 = f8
%即   C_a1a2 = ( A1_new + A2_new + A3_new + A4_new + A5_new + A6_new )*alpha_new + 3 - 2*( A1 + A2 + A3 );


g1 = diff(C_a1a2,a1);
g1 = collect(g1,[a1 a2]);
g1n = ( A_new(1,:) + A_new(2,:) + A_new(3,:) + A_new(4,:) + A_new(5,:) + A_new(6,:) ) * diff(alpha_new,a1) - 2*( A(1,:) + A(2,:) + A(3,:) )*diff(alpha,a1);
g1n = collect(g1n,[a1 a2]);
if ( simplify(g1 - g1n) == 0 )
    disp('12equal');  
end
%证得 g1 = ( A1_new + A2_new + A3_new + A4_new + A5_new + A6_new ) * diff(alpha_new,a1) - 2*( A1 + A2 + A3 )*diff(alpha,a1);
%同理 g2 = ( A1_new + A2_new + A3_new + A4_new + A5_new + A6_new ) * diff(alpha_new,a2) - 2*( A1 + A2 + A3 )*diff(alpha,a2);


%系数c11等真值
%分别对a1，a2求偏导
g1 = diff(C_a1a2,a1);
g2 = diff(C_a1a2,a2);
%使用coeffs函数提取系数，这里按照论文的方法将a2看作常数
coeff1 = coeffs(g1,a1);
coeff2 = coeffs(g2,a1);
c11 = coeff1(1,4);
c12 = coeff1(1,3);
c13 = coeff1(1,2);
c14 = coeff1(1,1);
c21 = coeff2(1,4);
c22 = coeff2(1,3);
c23 = coeff2(1,2);
c24 = coeff2(1,1);


Ans = ( A_new(1,:) + A_new(2,:) + A_new(3,:) + A_new(4,:) + A_new(5,:) + A_new(6,:) );
As = A(1,:) + A(2,:) + A(3,:);

coe11 = 4*Ans(1);
coe12 = 3*a2*Ans(2);
coe13 = 2*a2^2*Ans(3) + (-2)*2*As(1);
coe14 = a2^3*Ans(4) + (-2)*a2*As(2);
if ( coe11 - c11 == 0 )
    disp('13equal');  
end
if ( simplify(coe12 - c12) == 0 )
    disp('14equal');  
end
if ( simplify(coe13 - c13) == 0 )
    disp('15equal');  
end
if ( simplify(coe14 - c14) == 0 )
    disp('16equal');  
end
%证得   
% c11 = 4*Ans(1);
% c12 = 3*a2*Ans(2);
% c13 = 2*a2^2*Ans(3) + (-2)*2*As(1);
% c14 = a2^3*Ans(4) + (-2)*a2*As(2);
% Ans = A1_new + A2_new + A3_new + A4_new + A5_new + A6_new    1*5矩阵
% As = A(1,:) + A(2,:) + A(3,:)    1*3矩阵

coe21 = Ans(2);
coe22 = 2*a2*Ans(3);
coe23 = 3*a2^2*Ans(4) + (-2)*As(2);
coe24 = 4*a2^3*Ans(5) + (-2)*2*a2*As(3);
if ( simplify(coe21 - c21) == 0 )
    disp('17equal');  
end
if ( simplify(coe22 - c22) == 0 )
    disp('18equal');  
end
if ( simplify(coe23 - c23) == 0 )
    disp('19equal');  
end
if ( simplify(coe24 - c24) == 0 )
    disp('20equal');  
end
%证得   
% c21 = Ans(2);
% c22 = 2*a2*Ans(3);
% c23 = 3*a2^2*Ans(4) + (-2)*As(2);
% c24 = 4*a2^3*Ans(5) + (-2)*2*a2*As(3);
% Ans = A1_new + A2_new + A3_new + A4_new + A5_new + A6_new    1*5矩阵
% As = A(1,:) + A(2,:) + A(3,:)    1*3矩阵


%上面所证得的Ans、As为已知部分，可以由输入数据运算得出
%将Ans、As作为符号变量固定，求最终系数k9、k7等
syms Anss1 Anss2 Anss3 Anss4 Anss5 Ass1 Ass2 Ass3 real
Anss = [Anss1 Anss2 Anss3 Anss4 Anss5];
Ass = [Ass1 Ass2 Ass3];

coef11 = 4*Anss(1);
coef12 = 3*a2*Anss(2);
coef13 = 2*a2^2*Anss(3) + (-2)*2*Ass(1);
coef14 = a2^3*Anss(4) + (-2)*a2*Ass(2);

coef21 = Anss(2);
coef22 = 2*a2*Anss(3);
coef23 = 3*a2^2*Anss(4) + (-2)*Ass(2);
coef24 = 4*a2^3*Anss(5) + (-2)*2*a2*Ass(3);

mat1 = [coef11    0         0       coef21  0       0     ;
        coef12    coef11    0       coef22  coef21  0     ;
        coef13    coef12    coef11  coef23  coef22  coef21;
        coef14    coef13    coef12  coef24  coef23  coef22;
          0       coef14    coef13  0       coef24  coef23;
          0       0         coef14  0       0       coef24  ];

f_a = det(mat1);
[c,s] = coeffs(f_a,a2);
k9 = c(1);
k7 = c(2);
k5 = c(3);
k3 = c(4);
k1 = c(5);



