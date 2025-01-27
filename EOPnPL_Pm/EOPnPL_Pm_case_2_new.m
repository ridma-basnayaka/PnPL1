function [solution_r] = case_2_new( v1 , v2 )
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

coeffs_a2 = [k9 0 k7 0 k5 0 k3 0 k1];
solution_a2 = roots(coeffs_a2);
a2 = solution_a2( imag(solution_a2) == 0 );%取实根

%construct g1 expression
coe11 = 4*Ans(1);
coe12 = 3*a2*Ans(2);
coe13 = 2*a2.^2*Ans(3) + (-2)*2*As(1);
coe14 = a2.^3*Ans(4) + (-2)*a2*As(2);

n = size(a2,1);    %a2实数解个数
coeffs_a1 = [repmat(coe11,n,1) coe12 coe13 coe14];
a1a2 = [];
for i = 1 : n
   solution_a1 = roots(coeffs_a1(i,:)); %求解a1
   a1 = solution_a1( imag(solution_a1) == 0 );%取实根
   m = size(a1,1);   %k1实数解个数 
   solution_a1a2 = [a1 repmat(a2(i),m,1)];
   a1a2 = [a1a2 ; solution_a1a2 ];
end

solution_r = [v1 v2] * a1a2';
end


% syms a1 a2     
%求上面矩阵的行列式，得到函数f_a
% f_a = k9*a2^8 + k7*a2^6 + k5*a2^4 + k3*a2^2 + k1;
% solution_a2 = vpasolve(f_a,a2,[-inf inf]);

% g1 = coe11*a1^3 + coe12*a1^2 + coe13*a1 + coe14;

%substitute a2 into g1 to calculate a1
% generate a matrix where the right column contains solutions of a1 while
% corresponding solutions of a2 are on the left
% row = 1;
% for i = 1 : n
%     solution_g1 = subs(g1,a2,solution_a2(i));
%     solution = vpasolve(solution_g1,a1,[-inf inf]);
%     n1 = size(solution,1);
%     n_a2 = repmat(solution_a2(i),n1,1);
%     solution_a1a2(row : row + n1 - 1 , 1) = solution;
%     solution_a1a2(row : row + n1 - 1 , 2) = n_a2;
%     row = row + n1;
% end
% solution_r = [v1 v2] * solution_a1a2';  




