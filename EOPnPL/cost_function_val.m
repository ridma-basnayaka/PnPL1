function [cost_val] = cost_function_val(s,K)

s1 = s(1); s2 = s(2); s3 = s(3);

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
cost_val = dot(hr,ktkreshape);

   
end


