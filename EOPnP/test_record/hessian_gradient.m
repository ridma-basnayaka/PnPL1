clear
clc
syms s1 s2 s3 real; 
s = [s1;s2;s3];
s_x = [ 0 -s3 s2;
       s3  0 -s1;
       -s2 s1  0];
R_bar = (1-s'*s).*eye(3) + 2.*s_x + 2.*(s*s');
r_bar = reshape(R_bar',9,1);

ktk = sym('ktk',[9,9],'real');
Q = sym('Q',[9,9],'real');
cost_f = r_bar' * Q * r_bar;
H = hessian(cost_f,s);
% H = simplify(H)
G = gradient(cost_f,s);


hr1 = r_bar .* r_bar(1);    hr2 = r_bar .* r_bar(2);    hr3 = r_bar .* r_bar(3);    hr4 = r_bar .* r_bar(4);    hr5 = r_bar .* r_bar(5);
hr6 = r_bar .* r_bar(6);    hr7 = r_bar .* r_bar(7);    hr8 = r_bar .* r_bar(8);    hr9 = r_bar .* r_bar(9);
cf = dot(hr1,ktk(:,1)) + dot(hr2,ktk(:,2)) + dot(hr3,ktk(:,3)) + dot(hr4,ktk(:,4)) + dot(hr5,ktk(:,5)) + dot(hr6,ktk(:,6)) + ...
    dot(hr7,ktk(:,7)) + dot(hr8,ktk(:,8)) + dot(hr9,ktk(:,9));
if(simplify(cf - cost_f) == 0)
   disp('equal'); 
end

ktkreshape = reshape(ktk,1,81);
hr = [hr1;hr2;hr3;hr4;hr5;hr6;hr7;hr8;hr9];
cost = dot(hr,ktkreshape);
if(simplify(cost - cost_f) == 0)
   disp('equal'); 
end

hrdd11 = diff(hr,s1);
hrdd11 = diff(hrdd11,s1);
hrdd12 = diff(hr,s1);
hrdd12 = diff(hrdd12,s2);
hrdd13 = diff(hr,s1);
hrdd13 = diff(hrdd13,s3);

hrdd21 = diff(hr,s2);
hrdd21 = diff(hrdd21,s1);
hrdd22 = diff(hr,s2);
hrdd22 = diff(hrdd22,s2);
hrdd23 = diff(hr,s2);
hrdd23 = diff(hrdd23,s3);

hrdd31 = diff(hr,s3);
hrdd31 = diff(hrdd31,s1);
hrdd32 = diff(hr,s3);
hrdd32 = diff(hrdd32,s2);
hrdd33 = diff(hr,s3);
hrdd33 = diff(hrdd33,s3);

syms('hess',[3,3]);
hess(1,1) = dot(hrdd11,ktkreshape);
hess(1,2) = dot(hrdd12,ktkreshape);
hess(1,3) = dot(hrdd13,ktkreshape);
hess(2,1) = dot(hrdd21,ktkreshape);
hess(2,2) = dot(hrdd22,ktkreshape);
hess(2,3) = dot(hrdd23,ktkreshape);
hess(3,1) = dot(hrdd31,ktkreshape);
hess(3,2) = dot(hrdd32,ktkreshape);
hess(3,3) = dot(hrdd33,ktkreshape);

simplify(H - hess)


gr1 = diff(hr,s1);
gr2 = diff(hr,s2);
gr3 = diff(hr,s3);
g(1,1) = dot(gr1,ktkreshape);
g(2,1) = dot(gr2,ktkreshape);
g(3,1) = dot(gr3,ktkreshape);

simplify(G - g)