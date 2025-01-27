function [R] = CGR_degenerate()
syms a b c d real;
sum = a^2 + b^2 + c^2 + d^2;
abcd = [a b c d]';
quaternion_fractional = (1/sum) * [a^2+b^2-c^2-d^2, 2*b*c-2*a*d, 2*b*d+2*a*c;
                                 2*b*c+2*a*d, a^2-b^2+c^2-d^2, 2*c*d-2*a*b;
                                 2*b*d-2*a*c, 2*c*d+2*a*b, a^2-b^2-c^2+d^2];


ss = [0 1 2 4]';    %四元数四个参数
R = vpa( subs(quaternion_fractional,abcd,ss) );

end

