syms a b c d real;
s = a^2 + b^2 + c^2 + d^2;
quaternion_fractional = (1/s) * [a^2+b^2-c^2-d^2, 2*b*c-2*a*d, 2*b*d+2*a*c;
                                 2*b*c+2*a*d, a^2-b^2+c^2-d^2, 2*c*d-2*a*b;
                                 2*b*d-2*a*c, 2*c*d+2*a*b, a^2-b^2-c^2+d^2];
f1 = quaternion_fractional(2,1)*quaternion_fractional(3,2)-quaternion_fractional(2,2)*quaternion_fractional(3,1);
simplify(f1 - quaternion_fractional(1,3))
%                              
% syms w x y z real;
% s1 = [w x y z];
% quater = [1-2*y^2-2*z^2, 2*x*y+2*w*z, 2*x*z-2*w*y;
%               2*x*y-2*w*z, 1-2*x^2-2*z^2, 2*y*z+2*w*x;
%               2*x*z+2*w*y, 2*y*z-2*w*x, 1-2*x^2-2*y^2];

