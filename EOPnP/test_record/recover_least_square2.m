syms s1 s2 s3 real;
syms A B C real
s = [s1;s2;s3];
fs1 = (B^2+C^2+A^2) * s1^2 + 4*C*s1 + C^2;
fs2 = (B^2+C^2+A^2) * s2^2 - 4*B*s2 + B^2;
fs3 = (B^2+C^2+A^2) * s3^2 + 4*A*s3 + A^2;
cost_f = fs1^2 + fs2^2 + fs3^2;
f1 = diff(cost_f,s1);
f2 = diff(cost_f,s2);
f3 = diff(cost_f,s3);
% sol = solve([f1;f2;f3],s)
[cc1,ss1] = coeffs(f1,s);
[cc2,ss2] = coeffs(f2,s);
[cc3,ss3] = coeffs(f3,s);


