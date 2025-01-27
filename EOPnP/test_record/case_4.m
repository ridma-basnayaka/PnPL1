function [solution_r] = case_4(v1,v2,v3,v4)

%分解奇异向量
v11 = v1(1:3,:);    
v12 = v1(4:6,:);
v13 = v1(7:9,:);    %将奇异向量每三行分开，每个部分分别对应r1,r2,r3

v21 = v2(1:3,:);
v22 = v2(4:6,:);
v23 = v2(7:9,:);

v31 = v3(1:3,:);
v32 = v3(4:6,:);
v33 = v3(7:9,:);

v41 = v4(1:3,:);
v42 = v4(4:6,:);
v43 = v4(7:9,:);

syms a1 real
syms k1 k2 k3 real;

%旋转矩阵约束，顺便消去a1,a2,a3,a4
r1 = v11 + k1*v21 + k2*v31 + k3*v41;
r2 = v12 + k1*v22 + k2*v32 + k3*v42;
r3 = v13 + k1*v23 + k2*v33 + k3*v43;

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
% f1 = matlabFunction(M)
 
N_k = (-1) * [f4_k3k2(1,2) f5_k3k2(1,2) f6_k3k2(1,2);
              f4_k2k3(1,2) f5_k2k3(1,2) f6_k2k3(1,2);
              f4_k3k2(1,1) f5_k3k2(1,1) f6_k3k2(1,1)]';
% f2 = matlabFunction(N_k)
          
%求函数T_k
Tk = pinv(M) * N_k; 

%求M_k，根据等式（33）
% T_k(1,:) * [k2;k3;1] = k2^2
% T_k(2,:) * [k2;k3;1] = k3^2
% T_k(3,:) * [k2;k3;1] = k2 * k3
K = Tk * [k2;k3;1];
k2_2 = K(1,:);
k3_2 = K(2,:);
k2_k3 = K(3,:);

M_k(1,1) = k2_k3*Tk(1,1) + k3_2*Tk(1,2) + k3*Tk(1,3) - k2_2*Tk(3,1) - k2_k3*Tk(3,2) -k2*Tk(3,3);
M_k(2,1) = k2_k3*Tk(3,1) + k3_2*Tk(3,2) + k3*Tk(3,3) - k2_2*Tk(2,1) - k2_k3*Tk(2,2) -k2*Tk(2,3);
M_k(3,1) = k2_2*Tk(3,1)^2 + k2_k3*( 2*Tk(3,1)*Tk(3,2) ) + k3_2*Tk(3,2)^2 + 2*k2*Tk(3,1)*Tk(3,3) + 2*k3*Tk(3,2)*Tk(3,3) + Tk(3,3)^2 - ( k2_2*Tk(1,1)*Tk(2,1) + k2_k3*( Tk(1,1)*Tk(1,2) + Tk(1,2)*Tk(2,1) ) + k3_2*Tk(1,2)*Tk(2,2) + k2*( Tk(1,1)*Tk(2,3) + Tk(2,1)*Tk(2,3) )  + k3*( Tk(1,2)*Tk(2,3) + Tk(2,2)*Tk(1,3) ) + Tk(1,3)*Tk(2,3) );
%            ( k2_2*Tk(1,1)*Tk(2,1) + k2_k3*( Tk(1,1)*Tk(1,2) + Tk(1,2)*Tk(2,1) ) + k3_2*Tk(1,2)*Tk(2,2) + k2*( Tk(1,1)*Tk(2,3) + Tk(2,1)*Tk(2,3) )  + k3*( Tk(1,2)*Tk(2,3) + Tk(2,2)*Tk(1,3) ) + Tk(1,3)*Tk(2,3) ) 

%提取系数矩阵
[s1,~] = coeffs(M_k(1),[k2 k3]);
[s2,~] = coeffs(M_k(2),[k2 k3]);
[s3,~] = coeffs(M_k(3),[k2 k3]);
M_k = [s1;s2;s3];

% [cmk,smk] = coeffs(det(M_k));

%求行列式=0，求解k1
solution_k1 = vpasolve(det(M_k),k1,[-inf inf]);


%substitute k1 into M_k to get k2 and k3
%each column represents a solution(k1 k2 k3)
n_k1 = size(solution_k1,1);
for i = 1 : n_k1
    M_k = subs( M_k , k1 , solution_k1(i));
    M = M_k * [k2;k3;1];
    [solution_k2,solution_k3] = vpasolve( M(1:2) , [k2 k3] , [-inf inf] );
    solution_k123(:,i) = [solution_k1(i);solution_k2;solution_k3];
end

%use the cost function method introduced in the article to calculate a1
%different solutions of [k1 k2 k3] generate different r , cost and a1
% r = a1 .* ( [v1 v2 v3 v4] * [1;k2;k3;k4] )
r = a1 .* ( [v1 v2 v3 v4] * [ones(1,size(solution_k123,2));solution_k123] );
r1 = r(1:3,:);
r2 = r(4:6,:);
r3 = r(7:9,:);
cost_a1 = ( dot( r1 , r1) - 1 ).^2 + ( dot( r2 , r2 ) - 1 ).^2 + ( dot( r3 , r3 ) - 1 ).^2;
cost_a1 = diff( cost_a1 , a1 );

for i = 1 : n_k1
%     cost_a1 = ( dot( r1(:,i) , r1(:,i) ) - 1 )^2 + ( dot( r2(:,i) , r2(:,i) ) - 1 )^2 + ( dot( r3(:,i) , r3(:,i) ) - 1 )^2;
%     cost_a1 = diff( cost_a1 , a1 );
    v_a1 = vpasolve( cost_a1(i) , a1 , [-inf inf]); 
    
    %eliminate zero solution
    v_a1( find(v_a1==0) ) = [];
    
    solution_a1(:,i) = v_a1;    
end



n2 = size(solution_a1,1);
n3 = size(solution_k123,1);
n4 = size(solution_k123,2);
solution_a234 = zeros(n2*n3,n4);
if n2 == 1
   solution_a234 =  solution_a1 .* solution_k123;
else
    for i = 1 : n2
        for j = 1 : n3
           solution_a234(j + (i-1)*n3 , : ) = solution_a1(i,:) .*  solution_k123(j,:);
        end
    end
    solution_a234 = reshape(solution_a234, n3 , n2*n4 );
    solution_a1 = reshape(solution_a1, 1 , n2*n4 ); 
end



%each column is solution of [a1;a2;a3;a4]
solution_a1234 = [solution_a1;solution_a234];

%compute r result
solution_r = [v1 v2 v3 v4] * solution_a1234;

end

