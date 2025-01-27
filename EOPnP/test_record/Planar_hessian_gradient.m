clear
clc
syms s1 s2 s3 real; 
s = [s1;s2;s3];
s_x = [ 0 -s3 s2;
       s3  0 -s1;
       -s2 s1  0];
R_bar = (1-s'*s).*eye(3) + 2.*s_x + 2.*(s*s');
r_bar = reshape(R_bar',9,1);
r_bar = [r_bar(1);r_bar(2);r_bar(4);r_bar(5);r_bar(7);r_bar(8)];

Q = sym('Q',[6,6],'real');
cost_f = r_bar' * Q * r_bar;
H = hessian(cost_f,s);
% H = simplify(H)
G = gradient(cost_f,s);
% G = simplify(G)