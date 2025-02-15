syms k1 k2 k3 real
f(1,1) = 3*k1^2 + 4*k2^2 + 3*k3^2 + 5*k1*k2 + 2*k1*k3 + k2*k3 + k1 + k2 + 2*k3 - 14.42;
f(2,1) = 6*k1^2 + 4*k2^2 + 3*k3^2 + 8*k1*k2 + 25*k1*k3 + 7*k2*k3 + k1 + k2 + 2*k3 - 11.27;
f(3,1) = 4*k1^2 + 4*k2^2 + 3*k3^2 + 5*k1*k2 + 2*k1*k3 + 6*k2*k3 + k1 + k2 + 2*k3 - 12.42;
% k1 =0.5
% k2 = 1.5
% k3 = -0.3



f4_k3 = coeffs(f(1,1),k3);
f4_k3k2 = coeffs(f4_k3(1,1),k2);
f4_k2 = coeffs(f(1,1),k2);
f4_k2k3 = coeffs(f4_k2(1,1),k3);
f4_k2k3_2 = coeffs(f4_k2(1,2),k3);

f5_k3 = coeffs(f(2,1),k3);
f5_k3k2 = coeffs(f5_k3(1,1),k2);
f5_k2 = coeffs(f(2,1),k2);
f5_k2k3 = coeffs(f5_k2(1,1),k3);
f5_k2k3_2 = coeffs(f5_k2(1,2),k3);

f6_k3 = coeffs(f(3,1),k3);
f6_k3k2 = coeffs(f6_k3(1,1),k2);
f6_k2 = coeffs(f(3,1),k2);
f6_k2k3 = coeffs(f6_k2(1,1),k3);
f6_k2k3_2 = coeffs(f6_k2(1,2),k3);

M = [f4_k3k2(1,3) f4_k2k3(1,3) f4_k2k3_2(1,2);
     f5_k3k2(1,3) f5_k2k3(1,3) f5_k2k3_2(1,2);
     f6_k3k2(1,3) f6_k2k3(1,3) f6_k2k3_2(1,2)];
% f1 = matlabFunction(M)
 
N_k = (-1) * [f4_k3k2(1,2) f5_k3k2(1,2) f6_k3k2(1,2);
              f4_k2k3(1,2) f5_k2k3(1,2) f6_k2k3(1,2);
              f4_k3k2(1,1) f5_k3k2(1,1) f6_k3k2(1,1)]';

% M*[k2^2;k3^2;k2*k3]
% N_k*[k2;k3;1]


Tk = pinv(M) * N_k; 
K = Tk * [k2;k3;1];

k2_2 = K(1,:);
k3_2 = K(2,:);
k2_k3 = K(3,:);

M_k(1,1) = k2_k3*Tk(1,1) + k3_2*Tk(1,2) + k3*Tk(1,3) - k2_2*Tk(3,1) - k2_k3*Tk(3,2) -k2*Tk(3,3);
M_k(2,1) = k2_k3*Tk(3,1) + k3_2*Tk(3,2) + k3*Tk(3,3) - k2_2*Tk(2,1) - k2_k3*Tk(2,2) -k2*Tk(2,3);
M_k(3,1) = k2_2*Tk(3,1)^2 + k2_k3*( 2*Tk(3,1)*Tk(3,2) ) + k3_2*Tk(3,2)^2 + 2*k2*Tk(3,1)*Tk(3,3) + 2*k3*Tk(3,2)*Tk(3,3) + Tk(3,3)^2 - ( k2_2*Tk(1,1)*Tk(2,1) + k2_k3*( Tk(1,1)*Tk(1,2) + Tk(1,2)*Tk(2,1) ) + k3_2*Tk(1,2)*Tk(2,2) + k2*( Tk(1,1)*Tk(2,3) + Tk(2,1)*Tk(2,3) )  + k3*( Tk(1,2)*Tk(2,3) + Tk(2,2)*Tk(1,3) ) + Tk(1,3)*Tk(2,3) );
%            ( k2_2*Tk(1,1)*Tk(2,1) + k2_k3*( Tk(1,1)*Tk(1,2) + Tk(1,2)*Tk(2,1) ) + k3_2*Tk(1,2)*Tk(2,2) + k2*( Tk(1,1)*Tk(2,3) + Tk(2,1)*Tk(2,3) )  + k3*( Tk(1,2)*Tk(2,3) + Tk(2,2)*Tk(1,3) ) + Tk(1,3)*Tk(2,3) ) 

%提取系数矩阵
[s1,c1] = coeffs(M_k(1),[k2 k3]);
[s2,c2] = coeffs(M_k(2),[k2 k3]);
[s3,c3] = coeffs(M_k(3),[k2 k3]);
M_k = [s1;s2;s3];

%求行列式=0，求解k1
solution_k1 = vpasolve(det(M_k),k1,[-inf inf]);

