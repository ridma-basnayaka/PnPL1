function [s] = Planar_damped_newton(s0,K,W,V,Xw,Ximg,A)

n = size(s0,2);      %输入的s的个数
lambda = 10000*ones(1,n);
iter = 8;
Q = K'*K;
Q1_1 = Q(1,1);  Q1_2 = Q(1,2);  Q1_3 = Q(1,3);  Q1_4 = Q(1,4);  Q1_5 = Q(1,5);  Q1_6 = Q(1,6);  
Q2_1 = Q(2,1);  Q2_2 = Q(2,2);  Q2_3 = Q(2,3);  Q2_4 = Q(2,4);  Q2_5 = Q(2,5);  Q2_6 = Q(2,6);  
Q3_1 = Q(3,1);  Q3_2 = Q(3,2);  Q3_3 = Q(3,3);  Q3_4 = Q(3,4);  Q3_5 = Q(3,5);  Q3_6 = Q(3,6);   
Q4_1 = Q(4,1);  Q4_2 = Q(4,2);  Q4_3 = Q(4,3);  Q4_4 = Q(4,4);  Q4_5 = Q(4,5);  Q4_6 = Q(4,6);   
Q5_1 = Q(5,1);  Q5_2 = Q(5,2);  Q5_3 = Q(5,3);  Q5_4 = Q(5,4);  Q5_5 = Q(5,5);  Q5_6 = Q(5,6);   
Q6_1 = Q(6,1);  Q6_2 = Q(6,2);  Q6_3 = Q(6,3);  Q6_4 = Q(6,4);  Q6_5 = Q(6,5);  Q6_6 = Q(6,6);   

%coefficients for gradient
coeff1 = [ 4*Q1_1 - 4*Q1_4 - 4*Q4_1 + 4*Q4_4, 6*Q1_2 + 6*Q1_3 + 6*Q2_1 - 6*Q2_4 + 6*Q3_1 - 6*Q3_4 - 6*Q4_2 - 6*Q4_3, 6*Q1_5 - 6*Q4_5 + 6*Q5_1 - 6*Q5_4, 6*Q1_6 - 6*Q4_6 + 6*Q6_1 - 6*Q6_4, 4*Q1_4 - 4*Q1_1 + 8*Q2_2 + 8*Q2_3 + 8*Q3_2 + 8*Q3_3 + 4*Q4_1 - 4*Q4_4, 4*Q1_6 + 8*Q2_5 + 8*Q3_5 - 4*Q4_6 + 8*Q5_2 + 8*Q5_3 + 4*Q6_1 - 4*Q6_4, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 4*Q4_4 - 4*Q1_1 + 8*Q5_5, 4*Q1_3 - 4*Q1_2 - 4*Q2_1 + 4*Q2_4 + 4*Q3_1 - 4*Q3_4 + 4*Q4_2 - 4*Q4_3 + 8*Q5_6 + 8*Q6_5, 4*Q1_1 - 4*Q4_4 + 8*Q6_6, 2*Q2_4 - 2*Q1_3 - 2*Q2_1 - 2*Q1_2 - 2*Q3_1 + 2*Q3_4 + 2*Q4_2 + 2*Q4_3, 4*Q2_6 - 2*Q1_5 + 4*Q3_6 + 2*Q4_5 - 2*Q5_1 + 2*Q5_4 + 4*Q6_2 + 4*Q6_3, 2*Q4_6 - 4*Q2_5 - 4*Q3_5 - 2*Q1_6 - 4*Q5_2 - 4*Q5_3 - 2*Q6_1 + 2*Q6_4, 4*Q5_6 - 2*Q1_3 - 2*Q2_1 - 2*Q2_4 - 2*Q3_1 - 2*Q3_4 - 2*Q4_2 - 2*Q4_3 - 2*Q1_2 + 4*Q6_5, 8*Q3_3 - 8*Q2_2 - 8*Q5_5 + 8*Q6_6, 2*Q1_2 + 2*Q1_3 + 2*Q2_1 + 2*Q2_4 + 2*Q3_1 + 2*Q3_4 + 2*Q4_2 + 2*Q4_3 - 4*Q5_6 - 4*Q6_5, - 2*Q1_5 - 2*Q4_5 - 2*Q5_1 - 2*Q5_4, 4*Q3_5 - 4*Q2_5 - 2*Q1_6 - 2*Q4_6 - 4*Q5_2 + 4*Q5_3 - 2*Q6_1 - 2*Q6_4, 2*Q1_5 - 4*Q2_6 + 4*Q3_6 + 2*Q4_5 + 2*Q5_1 + 2*Q5_4 - 4*Q6_2 + 4*Q6_3, 2*Q1_6 + 2*Q4_6 + 2*Q6_1 + 2*Q6_4];
coeff2 = [ 2*Q1_2 + 2*Q1_3 + 2*Q2_1 - 2*Q2_4 + 2*Q3_1 - 2*Q3_4 - 2*Q4_2 - 2*Q4_3, 4*Q1_4 - 4*Q1_1 + 8*Q2_2 + 8*Q2_3 + 8*Q3_2 + 8*Q3_3 + 4*Q4_1 - 4*Q4_4, 2*Q1_6 + 4*Q2_5 + 4*Q3_5 - 2*Q4_6 + 4*Q5_2 + 4*Q5_3 + 2*Q6_1 - 2*Q6_4, 4*Q2_6 - 2*Q1_5 + 4*Q3_6 + 2*Q4_5 - 2*Q5_1 + 2*Q5_4 + 4*Q6_2 + 4*Q6_3, 6*Q2_4 - 6*Q1_3 - 6*Q2_1 - 6*Q1_2 - 6*Q3_1 + 6*Q3_4 + 6*Q4_2 + 6*Q4_3, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 4*Q4_6 - 8*Q2_5 - 8*Q3_5 - 4*Q1_6 - 8*Q5_2 - 8*Q5_3 - 4*Q6_1 + 4*Q6_4, 4*Q5_6 - 2*Q1_3 - 2*Q2_1 - 2*Q2_4 - 2*Q3_1 - 2*Q3_4 - 2*Q4_2 - 2*Q4_3 - 2*Q1_2 + 4*Q6_5, 8*Q3_3 - 8*Q2_2 - 8*Q5_5 + 8*Q6_6, 2*Q1_2 + 2*Q1_3 + 2*Q2_1 + 2*Q2_4 + 2*Q3_1 + 2*Q3_4 + 2*Q4_2 + 2*Q4_3 - 4*Q5_6 - 4*Q6_5, 4*Q1_1 - 4*Q1_4 - 4*Q4_1 + 4*Q4_4, 6*Q4_6 - 6*Q1_6 - 6*Q6_1 + 6*Q6_4, 6*Q1_5 - 6*Q4_5 + 6*Q5_1 - 6*Q5_4, 4*Q1_1 - 4*Q4_4 + 8*Q6_6, 4*Q1_2 - 4*Q1_3 + 4*Q2_1 - 4*Q2_4 - 4*Q3_1 + 4*Q3_4 - 4*Q4_2 + 4*Q4_3 - 8*Q5_6 - 8*Q6_5, 4*Q4_4 - 4*Q1_1 + 8*Q5_5, - 2*Q1_6 - 2*Q4_6 - 2*Q6_1 - 2*Q6_4, 2*Q1_5 - 4*Q2_6 + 4*Q3_6 + 2*Q4_5 + 2*Q5_1 + 2*Q5_4 - 4*Q6_2 + 4*Q6_3, 2*Q1_6 + 4*Q2_5 - 4*Q3_5 + 2*Q4_6 + 4*Q5_2 - 4*Q5_3 + 2*Q6_1 + 2*Q6_4, - 2*Q1_5 - 2*Q4_5 - 2*Q5_1 - 2*Q5_4];
coeff3 = [ 2*Q1_5 - 2*Q4_5 + 2*Q5_1 - 2*Q5_4, 2*Q1_6 + 4*Q2_5 + 4*Q3_5 - 2*Q4_6 + 4*Q5_2 + 4*Q5_3 + 2*Q6_1 - 2*Q6_4, 4*Q4_4 - 4*Q1_1 + 8*Q5_5, 2*Q1_3 - 2*Q1_2 - 2*Q2_1 + 2*Q2_4 + 2*Q3_1 - 2*Q3_4 + 2*Q4_2 - 2*Q4_3 + 4*Q5_6 + 4*Q6_5, 4*Q2_6 - 2*Q1_5 + 4*Q3_6 + 2*Q4_5 - 2*Q5_1 + 2*Q5_4 + 4*Q6_2 + 4*Q6_3, 8*Q5_6 - 4*Q1_3 - 4*Q2_1 - 4*Q2_4 - 4*Q3_1 - 4*Q3_4 - 4*Q4_2 - 4*Q4_3 - 4*Q1_2 + 8*Q6_5, 8*Q3_3 - 8*Q2_2 - 8*Q5_5 + 8*Q6_6, - 6*Q1_5 - 6*Q4_5 - 6*Q5_1 - 6*Q5_4, 8*Q3_5 - 8*Q2_5 - 4*Q1_6 - 4*Q4_6 - 8*Q5_2 + 8*Q5_3 - 4*Q6_1 - 4*Q6_4, 2*Q1_5 - 4*Q2_6 + 4*Q3_6 + 2*Q4_5 + 2*Q5_1 + 2*Q5_4 - 4*Q6_2 + 4*Q6_3, 2*Q4_6 - 2*Q1_6 - 2*Q6_1 + 2*Q6_4, 4*Q1_1 - 4*Q4_4 + 8*Q6_6, 2*Q1_2 - 2*Q1_3 + 2*Q2_1 - 2*Q2_4 - 2*Q3_1 + 2*Q3_4 - 2*Q4_2 + 2*Q4_3 - 4*Q5_6 - 4*Q6_5, - 6*Q1_6 - 6*Q4_6 - 6*Q6_1 - 6*Q6_4, 4*Q1_5 - 8*Q2_6 + 8*Q3_6 + 4*Q4_5 + 4*Q5_1 + 4*Q5_4 - 8*Q6_2 + 8*Q6_3, 2*Q1_6 + 4*Q2_5 - 4*Q3_5 + 2*Q4_6 + 4*Q5_2 - 4*Q5_3 + 2*Q6_1 + 2*Q6_4, 4*Q1_1 + 4*Q1_4 + 4*Q4_1 + 4*Q4_4, 6*Q1_2 - 6*Q1_3 + 6*Q2_1 + 6*Q2_4 - 6*Q3_1 - 6*Q3_4 + 6*Q4_2 - 6*Q4_3, 8*Q2_2 - 4*Q1_4 - 4*Q1_1 - 8*Q2_3 - 8*Q3_2 + 8*Q3_3 - 4*Q4_1 - 4*Q4_4, 2*Q1_3 - 2*Q1_2 - 2*Q2_1 - 2*Q2_4 + 2*Q3_1 + 2*Q3_4 - 2*Q4_2 + 2*Q4_3];

%coefficients for hessian
h11 = [ 12*Q1_1 - 12*Q1_4 - 12*Q4_1 + 12*Q4_4, 12*Q1_2 + 12*Q1_3 + 12*Q2_1 - 12*Q2_4 + 12*Q3_1 - 12*Q3_4 - 12*Q4_2 - 12*Q4_3, 12*Q1_5 - 12*Q4_5 + 12*Q5_1 - 12*Q5_4, 12*Q1_6 - 12*Q4_6 + 12*Q6_1 - 12*Q6_4, 4*Q1_4 - 4*Q1_1 + 8*Q2_2 + 8*Q2_3 + 8*Q3_2 + 8*Q3_3 + 4*Q4_1 - 4*Q4_4, 4*Q1_6 + 8*Q2_5 + 8*Q3_5 - 4*Q4_6 + 8*Q5_2 + 8*Q5_3 + 4*Q6_1 - 4*Q6_4, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 4*Q4_4 - 4*Q1_1 + 8*Q5_5, 4*Q1_3 - 4*Q1_2 - 4*Q2_1 + 4*Q2_4 + 4*Q3_1 - 4*Q3_4 + 4*Q4_2 - 4*Q4_3 + 8*Q5_6 + 8*Q6_5, 4*Q1_1 - 4*Q4_4 + 8*Q6_6];
h12 = [ 6*Q1_2 + 6*Q1_3 + 6*Q2_1 - 6*Q2_4 + 6*Q3_1 - 6*Q3_4 - 6*Q4_2 - 6*Q4_3, 8*Q1_4 - 8*Q1_1 + 16*Q2_2 + 16*Q2_3 + 16*Q3_2 + 16*Q3_3 + 8*Q4_1 - 8*Q4_4, 4*Q1_6 + 8*Q2_5 + 8*Q3_5 - 4*Q4_6 + 8*Q5_2 + 8*Q5_3 + 4*Q6_1 - 4*Q6_4, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 6*Q2_4 - 6*Q1_3 - 6*Q2_1 - 6*Q1_2 - 6*Q3_1 + 6*Q3_4 + 6*Q4_2 + 6*Q4_3, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 4*Q4_6 - 8*Q2_5 - 8*Q3_5 - 4*Q1_6 - 8*Q5_2 - 8*Q5_3 - 4*Q6_1 + 4*Q6_4, 4*Q5_6 - 2*Q1_3 - 2*Q2_1 - 2*Q2_4 - 2*Q3_1 - 2*Q3_4 - 2*Q4_2 - 2*Q4_3 - 2*Q1_2 + 4*Q6_5, 8*Q3_3 - 8*Q2_2 - 8*Q5_5 + 8*Q6_6, 2*Q1_2 + 2*Q1_3 + 2*Q2_1 + 2*Q2_4 + 2*Q3_1 + 2*Q3_4 + 2*Q4_2 + 2*Q4_3 - 4*Q5_6 - 4*Q6_5];
h13 = [ 6*Q1_5 - 6*Q4_5 + 6*Q5_1 - 6*Q5_4, 4*Q1_6 + 8*Q2_5 + 8*Q3_5 - 4*Q4_6 + 8*Q5_2 + 8*Q5_3 + 4*Q6_1 - 4*Q6_4, 8*Q4_4 - 8*Q1_1 + 16*Q5_5, 4*Q1_3 - 4*Q1_2 - 4*Q2_1 + 4*Q2_4 + 4*Q3_1 - 4*Q3_4 + 4*Q4_2 - 4*Q4_3 + 8*Q5_6 + 8*Q6_5, 4*Q2_6 - 2*Q1_5 + 4*Q3_6 + 2*Q4_5 - 2*Q5_1 + 2*Q5_4 + 4*Q6_2 + 4*Q6_3, 8*Q5_6 - 4*Q1_3 - 4*Q2_1 - 4*Q2_4 - 4*Q3_1 - 4*Q3_4 - 4*Q4_2 - 4*Q4_3 - 4*Q1_2 + 8*Q6_5, 8*Q3_3 - 8*Q2_2 - 8*Q5_5 + 8*Q6_6, - 6*Q1_5 - 6*Q4_5 - 6*Q5_1 - 6*Q5_4, 8*Q3_5 - 8*Q2_5 - 4*Q1_6 - 4*Q4_6 - 8*Q5_2 + 8*Q5_3 - 4*Q6_1 - 4*Q6_4, 2*Q1_5 - 4*Q2_6 + 4*Q3_6 + 2*Q4_5 + 2*Q5_1 + 2*Q5_4 - 4*Q6_2 + 4*Q6_3];

h22 = [ 4*Q1_4 - 4*Q1_1 + 8*Q2_2 + 8*Q2_3 + 8*Q3_2 + 8*Q3_3 + 4*Q4_1 - 4*Q4_4, 12*Q2_4 - 12*Q1_3 - 12*Q2_1 - 12*Q1_2 - 12*Q3_1 + 12*Q3_4 + 12*Q4_2 + 12*Q4_3, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 4*Q4_6 - 8*Q2_5 - 8*Q3_5 - 4*Q1_6 - 8*Q5_2 - 8*Q5_3 - 4*Q6_1 + 4*Q6_4, 12*Q1_1 - 12*Q1_4 - 12*Q4_1 + 12*Q4_4, 12*Q4_6 - 12*Q1_6 - 12*Q6_1 + 12*Q6_4, 12*Q1_5 - 12*Q4_5 + 12*Q5_1 - 12*Q5_4, 4*Q1_1 - 4*Q4_4 + 8*Q6_6, 4*Q1_2 - 4*Q1_3 + 4*Q2_1 - 4*Q2_4 - 4*Q3_1 + 4*Q3_4 - 4*Q4_2 + 4*Q4_3 - 8*Q5_6 - 8*Q6_5, 4*Q4_4 - 4*Q1_1 + 8*Q5_5];
h23 = [ 2*Q1_6 + 4*Q2_5 + 4*Q3_5 - 2*Q4_6 + 4*Q5_2 + 4*Q5_3 + 2*Q6_1 - 2*Q6_4, 8*Q2_6 - 4*Q1_5 + 8*Q3_6 + 4*Q4_5 - 4*Q5_1 + 4*Q5_4 + 8*Q6_2 + 8*Q6_3, 8*Q5_6 - 4*Q1_3 - 4*Q2_1 - 4*Q2_4 - 4*Q3_1 - 4*Q3_4 - 4*Q4_2 - 4*Q4_3 - 4*Q1_2 + 8*Q6_5, 8*Q3_3 - 8*Q2_2 - 8*Q5_5 + 8*Q6_6, 6*Q4_6 - 6*Q1_6 - 6*Q6_1 + 6*Q6_4, 8*Q1_1 - 8*Q4_4 + 16*Q6_6, 4*Q1_2 - 4*Q1_3 + 4*Q2_1 - 4*Q2_4 - 4*Q3_1 + 4*Q3_4 - 4*Q4_2 + 4*Q4_3 - 8*Q5_6 - 8*Q6_5, - 6*Q1_6 - 6*Q4_6 - 6*Q6_1 - 6*Q6_4, 4*Q1_5 - 8*Q2_6 + 8*Q3_6 + 4*Q4_5 + 4*Q5_1 + 4*Q5_4 - 8*Q6_2 + 8*Q6_3, 2*Q1_6 + 4*Q2_5 - 4*Q3_5 + 2*Q4_6 + 4*Q5_2 - 4*Q5_3 + 2*Q6_1 + 2*Q6_4];

h33 = [ 4*Q4_4 - 4*Q1_1 + 8*Q5_5, 8*Q5_6 - 4*Q1_3 - 4*Q2_1 - 4*Q2_4 - 4*Q3_1 - 4*Q3_4 - 4*Q4_2 - 4*Q4_3 - 4*Q1_2 + 8*Q6_5, - 12*Q1_5 - 12*Q4_5 - 12*Q5_1 - 12*Q5_4, 8*Q3_5 - 8*Q2_5 - 4*Q1_6 - 4*Q4_6 - 8*Q5_2 + 8*Q5_3 - 4*Q6_1 - 4*Q6_4, 4*Q1_1 - 4*Q4_4 + 8*Q6_6, - 12*Q1_6 - 12*Q4_6 - 12*Q6_1 - 12*Q6_4, 4*Q1_5 - 8*Q2_6 + 8*Q3_6 + 4*Q4_5 + 4*Q5_1 + 4*Q5_4 - 8*Q6_2 + 8*Q6_3, 12*Q1_1 + 12*Q1_4 + 12*Q4_1 + 12*Q4_4, 12*Q1_2 - 12*Q1_3 + 12*Q2_1 + 12*Q2_4 - 12*Q3_1 - 12*Q3_4 + 12*Q4_2 - 12*Q4_3, 8*Q2_2 - 4*Q1_4 - 4*Q1_1 - 8*Q2_3 - 8*Q3_2 + 8*Q3_3 - 4*Q4_1 - 4*Q4_4];

% strat the iteration
s = s0;
for i = 1 : n
%     cost_pre(i) = Planar_cost_function_val(s(:,i),K);
    [Re,te] = Planar_s2Rt(s(:,i),W,V);
    [err_pre(i)] = reprojection_error_usingRT(Xw',Ximg',Re,te,A);
end
eye = [1 0 0; 0 1 0; 0 0 1];
for i = 1 : iter 
    s1 = s(1,:);  s2 = s(2,:);  s3 = s(3,:);
    vec = [ s1.^2; s1.*s2; s1.*s3; s1; s2.^2; s2.*s3; s2; s3.^2; s3; ones(1,n)];
    vec3 = [s1.^3; s1.^2.*s2; s1.^2.*s3; s1.^2; s1.*s2.^2; s1.*s2.*s3; s1.*s2; s1.*s3.^2; s1.*s3; s1; s2.^3; s2.^2.*s3; s2.^2; s2.*s3.^2; s2.*s3; s2; s3.^3; s3.^2; s3; ones(1,n)];
    %gradient
    G = [coeff1*vec3; coeff2*vec3; coeff3*vec3];
    
    %hessian matrix
%     H = [h11*vec h12*vec h13*vec; 
%      h12*vec h22*vec h23*vec; 
%      h13*vec h23*vec h33*vec ];
    Hc1 = [h11*vec;h12*vec;h13*vec];
    Hc2 = [h12*vec;h22*vec;h23*vec];
    Hc3 = [h13*vec;h23*vec;h33*vec];
    
    for j = 1 : n
        d = ( [Hc1(:,j) Hc2(:,j) Hc3(:,j)] + lambda(:,j) * eye )\G(:,j);
        s(:,j) = s(:,j) - d;
%         cost_pre(i) = Planar_cost_function_val(s(:,i),K);
        [Re,te] = Planar_s2Rt(s(:,j),W,V);
        [err(:,j)] = reprojection_error_usingRT(Xw',Ximg',Re,te,A);
    end
   
    ind = (abs(err) < abs(err_pre));
    lambda(ind) = lambda(ind)/10;
    lambda(~ind) = lambda(~ind)*10;
    err_pre = err;
end

end

%     if abs(err) < abs(err_pre)
%         lambda = lambda/10;
%     else
%         lambda = lambda*10;
%     end
%     err_pre = err;   