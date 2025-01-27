syms s1 s2 s3 real;
ss = [s1;s2;s3];
s_x = [ 0 -s3 s2;
       s3  0 -s1;
       -s2 s1  0];
% Re = sym('Re',[3,3],'real');
% f1 = (Re - Re') - 4.*s_x/(1+ss'*ss);
% f2 = trace(Re) - (3-ss'*ss)/(1+ss'*ss);
% f3 = [f1(1,2);f1(1,3);f1(2,3)];
% sol_s = solve(f3);
% 
% A = Re(1,2) - Re(2,1);
% B = Re(1,3) - Re(3,1);
% C = Re(2,3) - Re(3,2);
% fs2 = (B^2+C^2+A^2) * s2^2 - 4*B*s2 + B^2;
% sol_s2 = solve(fs2);

syms A B C
f3(1) = A + (4*s3)/(s1^2 + s2^2 + s3^2 + 1);
f3(2) = B - (4*s2)/(s1^2 + s2^2 + s3^2 + 1);
f3(3) = C + (4*s1)/(s1^2 + s2^2 + s3^2 + 1);
sol_s = solve(f3);
s11 = sol_s.s1;
s22 = sol_s.s2;
s33 = sol_s.s3;

fs2 = (B^2+C^2+A^2) * s2^2 - 4*B*s2 + B^2;
sol_s2 = solve(fs2);

% simplify(s22 - sol_s2)

fs1 = (B^2+C^2+A^2) * s1^2 + 4*C*s1 + C^2;
sol_s1 = solve(fs1);

% simplify(s11 - sol_s1)

fs3 = (B^2+C^2+A^2) * s3^2 + 4*A*s3 + A^2;
sol_s3 = solve(fs3);

% simplify(s33 - sol_s3)

Re = sym('Re',[3,3],'real');   
Ai = Re(1,2) - Re(2,1);
Bi = Re(1,3) - Re(3,1);
Ci = Re(2,3) - Re(3,2);
fs3 = subs(fs3,[A B C],[Ai Bi Ci]);
fs2 = subs(fs2,[A B C],[Ai Bi Ci]);
fs1 = subs(fs1,[A B C],[Ai Bi Ci]);
[cs1,ss1] = coeffs(fs1,s1);
[cs2,ss2] = coeffs(fs2,s2);
[cs3,ss3] = coeffs(fs3,s3);






