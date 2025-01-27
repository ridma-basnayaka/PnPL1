function [R,t,R_bar,T] = Planar_s2Rt(s,W,V)
%由s生成对应的R和t

s_x = [ 0 -s(3) s(2);
       s(3)    0 -s(1);
       -s(2) s(1)  0];
R_bar = (1-s'*s).*eye(3) + 2.*s_x + 2.*(s*s');
R = R_bar / (1+s'*s);
r_bar = reshape(R_bar',9,1);
r_bar = [r_bar(1);r_bar(2);r_bar(4);r_bar(5);r_bar(7);r_bar(8)];
T = pinv(W) * V * r_bar;
t = T ./ (1 + s'*s);

end

