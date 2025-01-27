function [solution_r] = case_2( v1 , v2 )
%case_1()
%当r由两个右奇异向量表示时

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
f_a = det(mat1) / a2;
solution_a2 = vpasolve(f_a,a2,[-inf inf]);
n = size(solution_a2,1);    %number of solutions

%substitute a2 into g1 to calculate a1
% g1 = subs(g1,a2,solution_a2);

% generate a matrix where the right column contains solutions of a1 while
% corresponding solutions of a2 are on the left
row = 1;
for i = 1 : n
    solution_g1 = subs(g1,a2,solution_a2(i));
    solution = vpasolve(solution_g1,a1,[-inf inf]);
    n1 = size(solution,1);
    n_a2 = repmat(solution_a2(i),n1,1);
    solution_a1a2(row : row + n1 - 1 , 1) = solution;
    solution_a1a2(row : row + n1 - 1 , 2) = n_a2;
    row = row + n1;
end
solution_r = [v1 v2] * solution_a1a2';  


end

