function [s] = damped_newton_Ximgn(s0,K,W,V,Xw,Ximgn)
% damped newton method used in the article
% s0 ----- the initial solution s0

lambda = 10000;
iter = 5;
% initialize parameters
syms s1 s2 s3 real; 
ss = [s1;s2;s3];
r_bar = ...
   [s1*s1 - s2*s2 - s3*s3 + 1;
            2*s1*s2 - 2*s3;
            2*s2 + 2*s1*s3;
            2*s3 + 2*s1*s2;
  - s1*s1 + s2*s2 - s3*s3 + 1;
            2*s2*s3 - 2*s1;
            2*s1*s3 - 2*s2;
            2*s1 + 2*s2*s3;
  - s1*s1 - s2*s2 + s3*s3 + 1];
hr1 = r_bar .* r_bar(1);    hr2 = r_bar .* r_bar(2);    hr3 = r_bar .* r_bar(3);    hr4 = r_bar .* r_bar(4);    hr5 = r_bar .* r_bar(5);
hr6 = r_bar .* r_bar(6);    hr7 = r_bar .* r_bar(7);    hr8 = r_bar .* r_bar(8);    hr9 = r_bar .* r_bar(9);
ktk = K'*K;

ktkreshape = reshape(ktk,1,81);
hr = [hr1;hr2;hr3;hr4;hr5;hr6;hr7;hr8;hr9];

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

syms('Hess',[3,3]);
Hess(1,1) = dot(hrdd11,ktkreshape);
Hess(1,2) = dot(hrdd12,ktkreshape);
Hess(1,3) = dot(hrdd13,ktkreshape);
Hess(2,1) = dot(hrdd21,ktkreshape);
Hess(2,2) = dot(hrdd22,ktkreshape);
Hess(2,3) = dot(hrdd23,ktkreshape);
Hess(3,1) = dot(hrdd31,ktkreshape);
Hess(3,2) = dot(hrdd32,ktkreshape);
Hess(3,3) = dot(hrdd33,ktkreshape);
syms('Grad',[3,1]);
gr1 = diff(hr,s1);
gr2 = diff(hr,s2);
gr3 = diff(hr,s3);
Grad(1,1) = dot(gr1,ktkreshape);
Grad(2,1) = dot(gr2,ktkreshape);
Grad(3,1) = dot(gr3,ktkreshape);

% strat the iteration
s = s0;
[Re,te] = s2Rt(s,W,V);
[err_pre] = reprojection_error_Ximgn(Xw,Ximgn,Re,te);
eye = [1 0 0; 0 1 0; 0 0 1];
for i = 1 : iter 
    H = subs(Hess,ss,s);
    G = subs(Grad,ss,s);
    d = ( H + lambda * eye )\G;
    s = s - d;
    [Re,te] = s2Rt(s,W,V);
    [err] = reprojection_error_Ximgn(Xw,Ximgn,Re,te);

    if abs(err) < abs(err_pre)
        lambda = lambda/10;
    else
        lambda = lambda*10;
    end
    err_pre = err;          
end

end

