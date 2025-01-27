clear
clc
syms v11 v12 v13 v14 v15 v16 real;
syms v21 v22 v23 v24 v25 v26 real;
syms a1 a2 real;
v1 = [v11;v12;v13;v14;v15;v16];
v2 = [v21;v22;v23;v24;v25;v26];
r = a1*v1 + a2*v2;
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

unknowns = [a1 a2];
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
h8 = diff(cost_h,a2);
[ch7,sh7] = coeffs(h7,unknowns);
[ch8,sh8] = coeffs(h8,unknowns);
% disp(ch7);disp(ch8);
% disp(ch7(1:5));
% disp(ch7(6:10));
% disp(ch7(11:15));
% disp(ch7(16:20));
% 
% disp(ch8(1:5));
% disp(ch8(6:10));
% disp(ch8(11:15));
% disp(ch8(16:20));
% disp(sh7);disp(sh8);


m = sym('m',size(ch7),'real');
n = sym('n',size(ch8),'real');
% 隐变量方法 
c11 = m(1);    c12 = m(2)*a2;    c13 = m(3)*a2^2 + m(4);    c14 = m(5)*a2^3 + m(6)*a2;
c15 = m(7)*a2^4 + m(8)*a2^2 + m(9);    c16 = m(10)*a2^5 + m(11)*a2^3 + m(12)*a2;
c17 = m(13)*a2^6 + m(14)*a2^4 + m(15)*a2^2 + m(16);    c18 = m(17)*a2^7 + m(18)*a2^5 + m(19)*a2^3 + m(20)*a2;

c21 = n(1);    c22 = n(2)*a2;    c23 = n(3)*a2^2 + n(4);    c24 = n(5)*a2^3 + n(6)*a2;
c25 = n(7)*a2^4 + n(8)*a2^2 + n(9);    c26 = n(10)*a2^5 + n(11)*a2^3 + n(12)*a2;
c27 = n(13)*a2^6 + n(14)*a2^4 + n(15)*a2^2 + n(16);    c28 = n(17)*a2^7 + n(18)*a2^5 + n(19)*a2^3 + n(20)*a2;

% size of resultant: 2(x + 1) = 8 + x --> x = 6 
% --> size = (x+8)*(x+8) = 14*14
c1 = [c11;c12;c13;c14;c15;c16;c17;c18];
c2 = [c21;c22;c23;c24;c25;c26;c27;c28];
mat = sym(zeros(14,14));
for i = 1 : 7
   mat(i:i+7,i  ) = c1;
   mat(i:i+7,i+7) = c2;
end

% 直接求解输出的文本太长，此处分解 mat，并使用 polyeig（多项式特征值）求解
mat_0 = sym(zeros(14,14));
mat_1 = sym(zeros(14,14));
mat_2 = sym(zeros(14,14));
mat_3 = sym(zeros(14,14));
mat_4 = sym(zeros(14,14));
mat_5 = sym(zeros(14,14));
mat_6 = sym(zeros(14,14));
mat_7 = sym(zeros(14,14));

c1_0 = [m(1) 0 m(4) 0 m(9) 0 m(16) 0 ]';
c2_0 = [n(1) 0 n(4) 0 n(9) 0 n(16) 0 ]';

c1_1 = [0 m(2) 0 m(6) 0 m(12) 0 m(20)]';
c2_1 = [0 n(2) 0 n(6) 0 n(12) 0 n(20)]';

c1_2 = [0 0 m(3) 0 m(8) 0 m(15) 0]';
c2_2 = [0 0 n(3) 0 n(8) 0 n(15) 0]';

c1_3 = [0 0 0 m(5) 0 m(11) 0 m(19)]';
c2_3 = [0 0 0 n(5) 0 n(11) 0 n(19)]';

c1_4 = [0 0 0 0 m(7) 0 m(14) 0]';
c2_4 = [0 0 0 0 n(7) 0 n(14) 0]';

c1_5 = [0 0 0 0 0 m(10) 0 m(18)]';
c2_5 = [0 0 0 0 0 n(10) 0 n(18)]';

c1_6 = [0 0 0 0 0 0 m(13) 0]';
c2_6 = [0 0 0 0 0 0 n(13) 0]';

c1_7 = [0 0 0 0 0 0 0 m(17)]';
c2_7 = [0 0 0 0 0 0 0 n(17)]';


for i = 1 : 7
   mat_0(i:i+7,i  ) = c1_0;
   mat_0(i:i+7,i+7) = c2_0;
   mat_1(i:i+7,i  ) = c1_1;
   mat_1(i:i+7,i+7) = c2_1;
   mat_2(i:i+7,i  ) = c1_2;
   mat_2(i:i+7,i+7) = c2_2;
   mat_3(i:i+7,i  ) = c1_3;
   mat_3(i:i+7,i+7) = c2_3;
   mat_4(i:i+7,i  ) = c1_4;
   mat_4(i:i+7,i+7) = c2_4;
   mat_5(i:i+7,i  ) = c1_5;
   mat_5(i:i+7,i+7) = c2_5;
   mat_6(i:i+7,i  ) = c1_6;
   mat_6(i:i+7,i+7) = c2_6;
   mat_7(i:i+7,i  ) = c1_7;
   mat_7(i:i+7,i+7) = c2_7;
end



