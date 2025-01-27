clear
clc
syms v11 v12 v13 v14 v15 v16 real;
syms a1 real;
v1 = [v11;v12;v13;v14;v15;v16];
r = a1*v1;
c11 = r(1);
c21 = r(2);
c12 = r(3);
c22 = r(4);
c13 = r(5);
c23 = r(6);
c1 = [c11;c12;c13];
c2 = [c21;c22;c23];
c3 = cross(c1,c2);


% 利用行约束
c31 = c3(1);    c32 = c3(2);    c33 = c3(3);
r1 = [c11 c21 c31];
r2 = [c12 c22 c32];
r3 = [c13 c23 c33];

h1 = dot(r1,r1) - 1;
h2 = dot(r2,r2) - 1;
h3 = dot(r3,r3) - 1;
h4 = dot(r1,r2);
h5 = dot(r1,r3);
h6 = dot(r2,r3);

unknowns = a1;
[ch1,sh1] = coeffs(h1,unknowns);
[ch2,sh2] = coeffs(h2,unknowns);
[ch3,sh3] = coeffs(h3,unknowns);
[ch4,sh4] = coeffs(h4,unknowns);
[ch5,sh5] = coeffs(h5,unknowns);
[ch6,sh6] = coeffs(h6,unknowns);
% disp(ch1);disp(ch2);disp(ch3);disp(ch4);disp(ch5);disp(ch6);
% disp(sh1);disp(sh2);disp(sh3);disp(sh4);disp(sh5);disp(sh6);


cost_h = h1^2 + h2^2 + h3^2 + h4^2 + h5^2 + h6^2;
h7 = diff(cost_h,a1);
[ch7,sh7] = coeffs(h7,a1);
disp(ch7);

% 化简
sh7 = sh7./a1;
disp(sh7);





