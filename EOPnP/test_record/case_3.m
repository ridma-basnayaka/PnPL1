function [solution_r] = case_3(v1,v2,v3)

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

%设置系数未知数
syms a1 a2 a3 real; 
syms k1 k2 real;

%构建三个r的等式，并利用约束，顺便消去a1
r1 = (v11 + k1.*v21 + k2.*v31);
r2 = (v12 + k1.*v22 + k2.*v32);
r3 = (v13 + k1.*v23 + k2.*v33);

%矩阵B元素
g1 = dot(r1,r1) - dot(r2,r2);
g2 = dot(r1,r1) - dot(r3,r3);
f4 = dot(r1,r2);
f5 = dot(r1,r3);
f6 = dot(r2,r3);
% B = [g1;g2;f4;f5;f6];

%求代价函数C_k1k2
C_k1k2 = (g1)^2 + (g2)^2 + (f4)^2 + (f5)^2 + (f6)^2;



%分别对k1，k2求偏导
t1 = diff(C_k1k2,k1);
t2 = diff(C_k1k2,k2);


%使用coeffs函数提取系数，这里按照论文的方法将k2看作常数
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

%构造文中（22）矩阵
mat1 = [c11   0   0 c21   0   0;
        c12 c11   0 c22 c21   0;
        c13 c12 c11 c23 c22 c21;
        c14 c13 c12 c24 c23 c22;
          0 c14 c13   0 c24 c23
          0   0 c14   0   0 c24];
      
%求上面矩阵的行列式，得到函数f_k,
f_k = det(mat1);
solution_k2 = vpasolve(f_k,k2,[-inf inf]);    %添加[-inf inf]，是函数只取实数解




% generate a matrix where the right column contains solutions of k1 while
% corresponding solutions of k2 are on the left
n = size(solution_k2,1);    %number of solutions

%substitute k2 into t1 to calculate k1
t1 = subs(t1,k2,solution_k2);
row = 1;
for i = 1 : n
    solution = vpasolve(t1(i),k1,[-inf inf]);
    n1 = size(solution,1);
    n = repmat(solution_k2(i),n1,1);
    solution_k12(row : row + n1 - 1 , 1) = solution;
    solution_k12(row : row + n1 - 1 , 2) = n;
    row = row + n1;
end



% substitute above equations into f1 = 1 to compute a1
% f1 = dot(r1,r1) = 1
% r = a1 .* [v1 v2 v3] * [1;k1;k2]
r = a1 .* [v1 v2 v3] * [ones(size(solution_k12,1),1) solution_k12]';
r1 = r(1:3,:);
f1 = dot(r1,r1) - 1;    %every solution of k1 k2 generates a f1(every row) and then generate corresponding a1

for i = 1 : size(f1,2)
    solution_a1(:,i) = vpasolve(f1(i),a1,[-inf inf]);  
end


solution_k12 = solution_k12';
n2 = size(solution_a1,1);
n3 = size(solution_k12,1);
n4 = size(solution_k12,2);
for i = 1 : n2
    for j = 1 : n3
       solution_a23(j + (i-1)*n3 , : ) = solution_a1(i,:) .*  solution_k12(j,:);
    end
    
end

% % although the function reshape() reshape the matrix according to
% % column,but here we reshape them in the same manner so that we do not
% % destroy the correspondences of a1 a2 a3
solution_a23 = reshape(solution_a23, n3 , n2*n4 );
solution_a1 = reshape(solution_a1, 1 , n2*n4 );

%each column is solution of [a1;a2;a3]
solution_a123 = [solution_a1;solution_a23];
% solution_a123 = [reshape(solution_a1', 1 ,n2*n4 );
%                 solution_a23(1:n3,1:n4)];

%compute r result
solution_r = [v1 v2 v3] * solution_a123;




end

