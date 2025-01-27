syms Ans1 Ans2 Ans3 Ans4 Ans5 a1 a2 As1 As2 As3 real
c11 = 4*Ans1;
c12 = 3*a2*Ans2;
c13 = 2*a2^2 + 2*As1*(-2);
c14 = a2^3*Ans4 + a2*As2*(-2);

c21 = Ans2;
c22 = 2*a2*Ans3;
c23 = 3*a2^2*Ans4 + As2*(-2);
c24 = 4*a2^3*Ans5 + 2*a2*As3*(-2);

c1 = [c11 c12 c13 c14];
c2 = [c21 c22 c23 c24];

% c3 = c1 - (c14/c24)*c2;  %c3(4) = 0
% f1 = c3 * [a1^3 a1^2 a1 1]';
% [cf1,sf1] = coeffs(f1,a1);

% original polynomial equations system
f1 = c1 * [a1^3 a1^2 a1 1]';
f2 = c2 * [a1^3 a1^2 a1 1]';






