function [H,G,f] = EOPnP_compute_LM_parameter(si,K)

% construct the cost function first
syms s1 s2 s3 real;
s = [s1;s2;s3];
s_x = [ 0 -s3 s2;
       s3  0 -s1;
       -s2 s1  0];
R_bar = (1-s'*s).*eye(3) + 2.*s_x + 2.*(s*s');
r_bar = reshape(R_bar',9,1);
cost_f = r_bar' * (K' * K) * r_bar;

% compute gradient vector and hessian matrix
G = gradient(cost_f);
H = hessian(cost_f);

%compute H and G and cost function
G = subs(G,s,si);
H = subs(H,s,si);

% compute the error, since the ground truth of the cost function is 
% equal to zero,the error is equal to the result of the cost function
f = subs(cost_f,s,si);



end

