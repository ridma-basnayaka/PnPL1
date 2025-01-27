syms s1 s2 s3 real; 
re = sym('re',[9,1],'real');
s = [s1;s2;s3];
s_x = [ 0 -s3 s2;
       s3  0 -s1;
       -s2 s1  0];
R_bar = (1-s'*s).*eye(3) + 2.*s_x + 2.*(s*s');
R = R_bar / (1+s'*s);
r = reshape(R',9,1);
f = r - re;
cost_f = f(1)^2 + f(2)^2 + f(3)^2 + f(4)^2 + f(5)^2 + f(6)^2 + f(7)^2 + f(8)^2 + f(9)^2;
f1 = simplify(diff(cost_f,s1));
f2 = simplify(diff(cost_f,s2));
f3 = simplify(diff(cost_f,s3));

%等于0，因此可以消去系数 (1+s'*s)^2 
f1 = simplify(f1 * (1+s'*s)^2);
f2 = simplify(f2 * (1+s'*s)^2);
f3 = simplify(f3 * (1+s'*s)^2);
% solve([f1;f2;f3],s)

