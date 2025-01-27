function [s] = Polish_s_reduce_cost(s0,K)
% damped newton method used in the article
% s0 ----- the initial solution s0

lambda = 10000;
iter = 8;
Q = K'*K;
Q1_1 = Q(1,1);  Q1_2 = Q(1,2);  Q1_3 = Q(1,3);  Q1_4 = Q(1,4);  Q1_5 = Q(1,5);  Q1_6 = Q(1,6);  Q1_7 = Q(1,7);  Q1_8 = Q(1,8);  Q1_9 = Q(1,9);  
Q2_1 = Q(2,1);  Q2_2 = Q(2,2);  Q2_3 = Q(2,3);  Q2_4 = Q(2,4);  Q2_5 = Q(2,5);  Q2_6 = Q(2,6);  Q2_7 = Q(2,7);  Q2_8 = Q(2,8);  Q2_9 = Q(2,9);  
Q3_1 = Q(3,1);  Q3_2 = Q(3,2);  Q3_3 = Q(3,3);  Q3_4 = Q(3,4);  Q3_5 = Q(3,5);  Q3_6 = Q(3,6);  Q3_7 = Q(3,7);  Q3_8 = Q(3,8);  Q3_9 = Q(3,9);  
Q4_1 = Q(4,1);  Q4_2 = Q(4,2);  Q4_3 = Q(4,3);  Q4_4 = Q(4,4);  Q4_5 = Q(4,5);  Q4_6 = Q(4,6);  Q4_7 = Q(4,7);  Q4_8 = Q(4,8);  Q4_9 = Q(4,9);  
Q5_1 = Q(5,1);  Q5_2 = Q(5,2);  Q5_3 = Q(5,3);  Q5_4 = Q(5,4);  Q5_5 = Q(5,5);  Q5_6 = Q(5,6);  Q5_7 = Q(5,7);  Q5_8 = Q(5,8);  Q5_9 = Q(5,9); 
Q6_1 = Q(6,1);  Q6_2 = Q(6,2);  Q6_3 = Q(6,3);  Q6_4 = Q(6,4);  Q6_5 = Q(6,5);  Q6_6 = Q(6,6);  Q6_7 = Q(6,7);  Q6_8 = Q(6,8);  Q6_9 = Q(6,9); 
Q7_1 = Q(7,1);  Q7_2 = Q(7,2);  Q7_3 = Q(7,3);  Q7_4 = Q(7,4);  Q7_5 = Q(7,5);  Q7_6 = Q(7,6);  Q7_7 = Q(7,7);  Q7_8 = Q(7,8);  Q7_9 = Q(7,9); 
Q8_1 = Q(8,1);  Q8_2 = Q(8,2);  Q8_3 = Q(8,3);  Q8_4 = Q(8,4);  Q8_5 = Q(8,5);  Q8_6 = Q(8,6);  Q8_7 = Q(8,7);  Q8_8 = Q(8,8);  Q8_9 = Q(8,9); 
Q9_1 = Q(9,1);  Q9_2 = Q(1,2);  Q9_3 = Q(9,3);  Q9_4 = Q(9,4);  Q9_5 = Q(9,5);  Q9_6 = Q(9,6);  Q9_7 = Q(9,7);  Q9_8 = Q(9,8);  Q9_9 = Q(9,9); 

%coefficients for gradient
coeff1 = [ 4*Q1_1 - 4*Q1_5 - 4*Q1_9 - 4*Q5_1 + 4*Q5_5 + 4*Q5_9 - 4*Q9_1 + 4*Q9_5 + 4*Q9_9, 6*Q1_2 + 6*Q1_4 + 6*Q2_1 - 6*Q2_5 - 6*Q2_9 + 6*Q4_1 - 6*Q4_5 - 6*Q4_9 - 6*Q5_2 - 6*Q5_4 - 6*Q9_2 - 6*Q9_4, 6*Q1_3 + 6*Q1_7 + 6*Q3_1 - 6*Q3_5 - 6*Q3_9 - 6*Q5_3 - 6*Q5_7 + 6*Q7_1 - 6*Q7_5 - 6*Q7_9 - 6*Q9_3 - 6*Q9_7, 6*Q1_8 - 6*Q1_6 + 6*Q5_6 - 6*Q5_8 - 6*Q6_1 + 6*Q6_5 + 6*Q6_9 + 6*Q8_1 - 6*Q8_5 - 6*Q8_9 + 6*Q9_6 - 6*Q9_8, 4*Q1_5 - 4*Q1_1 + 8*Q2_2 + 8*Q2_4 + 8*Q4_2 + 8*Q4_4 + 4*Q5_1 - 4*Q5_5 + 4*Q9_9, 4*Q1_6 + 4*Q1_8 + 8*Q2_3 + 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 - 4*Q5_6 - 4*Q5_8 + 4*Q6_1 - 4*Q6_5 - 4*Q6_9 + 8*Q7_2 + 8*Q7_4 + 4*Q8_1 - 4*Q8_5 - 4*Q8_9 - 4*Q9_6 - 4*Q9_8, 4*Q1_3 - 4*Q1_7 - 8*Q2_6 + 8*Q2_8 + 4*Q3_1 - 4*Q3_5 - 4*Q3_9 - 8*Q4_6 + 8*Q4_8 - 4*Q5_3 + 4*Q5_7 - 8*Q6_2 - 8*Q6_4 - 4*Q7_1 + 4*Q7_5 + 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 + 4*Q9_7, 4*Q1_9 - 4*Q1_1 + 8*Q3_3 + 8*Q3_7 + 4*Q5_5 + 8*Q7_3 + 8*Q7_7 + 4*Q9_1 - 4*Q9_9, 4*Q1_4 - 4*Q1_2 - 4*Q2_1 + 4*Q2_5 + 4*Q2_9 - 8*Q3_6 + 8*Q3_8 + 4*Q4_1 - 4*Q4_5 - 4*Q4_9 + 4*Q5_2 - 4*Q5_4 - 8*Q6_3 - 8*Q6_7 - 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 - 4*Q9_4, 4*Q1_1 - 4*Q5_5 - 4*Q5_9 + 8*Q6_6 - 8*Q6_8 - 8*Q8_6 + 8*Q8_8 - 4*Q9_5 - 4*Q9_9, 2*Q2_5 - 2*Q1_4 - 2*Q2_1 - 2*Q1_2 - 2*Q2_9 - 2*Q4_1 + 2*Q4_5 - 2*Q4_9 + 2*Q5_2 + 2*Q5_4 - 2*Q9_2 - 2*Q9_4, 4*Q2_6 - 2*Q1_7 - 2*Q1_3 + 4*Q2_8 - 2*Q3_1 + 2*Q3_5 - 2*Q3_9 + 4*Q4_6 + 4*Q4_8 + 2*Q5_3 + 2*Q5_7 + 4*Q6_2 + 4*Q6_4 - 2*Q7_1 + 2*Q7_5 - 2*Q7_9 + 4*Q8_2 + 4*Q8_4 - 2*Q9_3 - 2*Q9_7, 2*Q1_6 - 2*Q1_8 + 4*Q2_3 - 4*Q2_7 + 4*Q3_2 + 4*Q3_4 + 4*Q4_3 - 4*Q4_7 - 2*Q5_6 + 2*Q5_8 + 2*Q6_1 - 2*Q6_5 + 2*Q6_9 - 4*Q7_2 - 4*Q7_4 - 2*Q8_1 + 2*Q8_5 - 2*Q8_9 + 2*Q9_6 - 2*Q9_8, 2*Q2_9 - 2*Q1_4 - 2*Q2_1 - 2*Q2_5 - 2*Q1_2 + 4*Q3_6 + 4*Q3_8 - 2*Q4_1 - 2*Q4_5 + 2*Q4_9 - 2*Q5_2 - 2*Q5_4 + 4*Q6_3 + 4*Q6_7 + 4*Q7_6 + 4*Q7_8 + 4*Q8_3 + 4*Q8_7 + 2*Q9_2 + 2*Q9_4, 8*Q3_3 - 8*Q2_2 + 8*Q4_4 - 8*Q6_6 - 8*Q7_7 + 8*Q8_8, 2*Q1_2 + 2*Q1_4 + 2*Q2_1 + 2*Q2_5 + 2*Q2_9 - 4*Q3_6 + 4*Q3_8 + 2*Q4_1 + 2*Q4_5 + 2*Q4_9 + 2*Q5_2 + 2*Q5_4 - 4*Q6_3 + 4*Q6_7 + 4*Q7_6 - 4*Q7_8 + 4*Q8_3 - 4*Q8_7 + 2*Q9_2 + 2*Q9_4, 2*Q3_9 - 2*Q1_7 - 2*Q3_1 - 2*Q3_5 - 2*Q1_3 - 2*Q5_3 - 2*Q5_7 - 2*Q7_1 - 2*Q7_5 + 2*Q7_9 + 2*Q9_3 + 2*Q9_7, 2*Q1_6 - 2*Q1_8 - 4*Q2_3 - 4*Q2_7 - 4*Q3_2 + 4*Q3_4 + 4*Q4_3 + 4*Q4_7 + 2*Q5_6 - 2*Q5_8 + 2*Q6_1 + 2*Q6_5 - 2*Q6_9 - 4*Q7_2 + 4*Q7_4 - 2*Q8_1 - 2*Q8_5 + 2*Q8_9 - 2*Q9_6 + 2*Q9_8, 2*Q1_3 + 2*Q1_7 + 4*Q2_6 - 4*Q2_8 + 2*Q3_1 + 2*Q3_5 + 2*Q3_9 - 4*Q4_6 + 4*Q4_8 + 2*Q5_3 + 2*Q5_7 + 4*Q6_2 - 4*Q6_4 + 2*Q7_1 + 2*Q7_5 + 2*Q7_9 - 4*Q8_2 + 4*Q8_4 + 2*Q9_3 + 2*Q9_7, 2*Q1_8 - 2*Q1_6 - 2*Q5_6 + 2*Q5_8 - 2*Q6_1 - 2*Q6_5 - 2*Q6_9 + 2*Q8_1 + 2*Q8_5 + 2*Q8_9 - 2*Q9_6 + 2*Q9_8];
coeff2 = [ 2*Q1_2 + 2*Q1_4 + 2*Q2_1 - 2*Q2_5 - 2*Q2_9 + 2*Q4_1 - 2*Q4_5 - 2*Q4_9 - 2*Q5_2 - 2*Q5_4 - 2*Q9_2 - 2*Q9_4, 4*Q1_5 - 4*Q1_1 + 8*Q2_2 + 8*Q2_4 + 8*Q4_2 + 8*Q4_4 + 4*Q5_1 - 4*Q5_5 + 4*Q9_9, 2*Q1_6 + 2*Q1_8 + 4*Q2_3 + 4*Q2_7 + 4*Q3_2 + 4*Q3_4 + 4*Q4_3 + 4*Q4_7 - 2*Q5_6 - 2*Q5_8 + 2*Q6_1 - 2*Q6_5 - 2*Q6_9 + 4*Q7_2 + 4*Q7_4 + 2*Q8_1 - 2*Q8_5 - 2*Q8_9 - 2*Q9_6 - 2*Q9_8, 2*Q1_3 - 2*Q1_7 - 4*Q2_6 + 4*Q2_8 + 2*Q3_1 - 2*Q3_5 - 2*Q3_9 - 4*Q4_6 + 4*Q4_8 - 2*Q5_3 + 2*Q5_7 - 4*Q6_2 - 4*Q6_4 - 2*Q7_1 + 2*Q7_5 + 2*Q7_9 + 4*Q8_2 + 4*Q8_4 - 2*Q9_3 + 2*Q9_7, 6*Q2_5 - 6*Q1_4 - 6*Q2_1 - 6*Q1_2 - 6*Q2_9 - 6*Q4_1 + 6*Q4_5 - 6*Q4_9 + 6*Q5_2 + 6*Q5_4 - 6*Q9_2 - 6*Q9_4, 8*Q2_6 - 4*Q1_7 - 4*Q1_3 + 8*Q2_8 - 4*Q3_1 + 4*Q3_5 - 4*Q3_9 + 8*Q4_6 + 8*Q4_8 + 4*Q5_3 + 4*Q5_7 + 8*Q6_2 + 8*Q6_4 - 4*Q7_1 + 4*Q7_5 - 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 - 4*Q9_7, 4*Q1_6 - 4*Q1_8 + 8*Q2_3 - 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 - 8*Q4_7 - 4*Q5_6 + 4*Q5_8 + 4*Q6_1 - 4*Q6_5 + 4*Q6_9 - 8*Q7_2 - 8*Q7_4 - 4*Q8_1 + 4*Q8_5 - 4*Q8_9 + 4*Q9_6 - 4*Q9_8, 2*Q2_9 - 2*Q1_4 - 2*Q2_1 - 2*Q2_5 - 2*Q1_2 + 4*Q3_6 + 4*Q3_8 - 2*Q4_1 - 2*Q4_5 + 2*Q4_9 - 2*Q5_2 - 2*Q5_4 + 4*Q6_3 + 4*Q6_7 + 4*Q7_6 + 4*Q7_8 + 4*Q8_3 + 4*Q8_7 + 2*Q9_2 + 2*Q9_4, 8*Q3_3 - 8*Q2_2 + 8*Q4_4 - 8*Q6_6 - 8*Q7_7 + 8*Q8_8, 2*Q1_2 + 2*Q1_4 + 2*Q2_1 + 2*Q2_5 + 2*Q2_9 - 4*Q3_6 + 4*Q3_8 + 2*Q4_1 + 2*Q4_5 + 2*Q4_9 + 2*Q5_2 + 2*Q5_4 - 4*Q6_3 + 4*Q6_7 + 4*Q7_6 - 4*Q7_8 + 4*Q8_3 - 4*Q8_7 + 2*Q9_2 + 2*Q9_4, 4*Q1_1 - 4*Q1_5 + 4*Q1_9 - 4*Q5_1 + 4*Q5_5 - 4*Q5_9 + 4*Q9_1 - 4*Q9_5 + 4*Q9_9, 6*Q5_6 - 6*Q1_8 - 6*Q1_6 + 6*Q5_8 - 6*Q6_1 + 6*Q6_5 - 6*Q6_9 - 6*Q8_1 + 6*Q8_5 - 6*Q8_9 - 6*Q9_6 - 6*Q9_8, 6*Q1_7 - 6*Q1_3 - 6*Q3_1 + 6*Q3_5 - 6*Q3_9 + 6*Q5_3 - 6*Q5_7 + 6*Q7_1 - 6*Q7_5 + 6*Q7_9 - 6*Q9_3 + 6*Q9_7, 4*Q1_1 - 4*Q5_5 + 4*Q5_9 + 8*Q6_6 + 8*Q6_8 + 8*Q8_6 + 8*Q8_8 + 4*Q9_5 - 4*Q9_9, 4*Q1_2 - 4*Q1_4 + 4*Q2_1 - 4*Q2_5 + 4*Q2_9 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 + 4*Q4_5 - 4*Q4_9 - 4*Q5_2 + 4*Q5_4 + 8*Q6_3 - 8*Q6_7 - 8*Q7_6 - 8*Q7_8 + 8*Q8_3 - 8*Q8_7 + 4*Q9_2 - 4*Q9_4, 8*Q3_3 - 4*Q1_9 - 4*Q1_1 - 8*Q3_7 + 4*Q5_5 - 8*Q7_3 + 8*Q7_7 - 4*Q9_1 - 4*Q9_9, 2*Q6_9 - 2*Q1_8 - 2*Q5_6 - 2*Q5_8 - 2*Q6_1 - 2*Q6_5 - 2*Q1_6 - 2*Q8_1 - 2*Q8_5 + 2*Q8_9 + 2*Q9_6 + 2*Q9_8, 2*Q1_7 - 2*Q1_3 - 4*Q2_6 - 4*Q2_8 - 2*Q3_1 - 2*Q3_5 + 2*Q3_9 + 4*Q4_6 + 4*Q4_8 - 2*Q5_3 + 2*Q5_7 - 4*Q6_2 + 4*Q6_4 + 2*Q7_1 + 2*Q7_5 - 2*Q7_9 - 4*Q8_2 + 4*Q8_4 + 2*Q9_3 - 2*Q9_7, 2*Q1_6 + 2*Q1_8 - 4*Q2_3 + 4*Q2_7 - 4*Q3_2 + 4*Q3_4 + 4*Q4_3 - 4*Q4_7 + 2*Q5_6 + 2*Q5_8 + 2*Q6_1 + 2*Q6_5 + 2*Q6_9 + 4*Q7_2 - 4*Q7_4 + 2*Q8_1 + 2*Q8_5 + 2*Q8_9 + 2*Q9_6 + 2*Q9_8, 2*Q1_3 - 2*Q1_7 + 2*Q3_1 + 2*Q3_5 + 2*Q3_9 + 2*Q5_3 - 2*Q5_7 - 2*Q7_1 - 2*Q7_5 - 2*Q7_9 + 2*Q9_3 - 2*Q9_7];
coeff3 = [ 2*Q1_3 + 2*Q1_7 + 2*Q3_1 - 2*Q3_5 - 2*Q3_9 - 2*Q5_3 - 2*Q5_7 + 2*Q7_1 - 2*Q7_5 - 2*Q7_9 - 2*Q9_3 - 2*Q9_7, 2*Q1_6 + 2*Q1_8 + 4*Q2_3 + 4*Q2_7 + 4*Q3_2 + 4*Q3_4 + 4*Q4_3 + 4*Q4_7 - 2*Q5_6 - 2*Q5_8 + 2*Q6_1 - 2*Q6_5 - 2*Q6_9 + 4*Q7_2 + 4*Q7_4 + 2*Q8_1 - 2*Q8_5 - 2*Q8_9 - 2*Q9_6 - 2*Q9_8, 4*Q1_9 - 4*Q1_1 + 8*Q3_3 + 8*Q3_7 + 4*Q5_5 + 8*Q7_3 + 8*Q7_7 + 4*Q9_1 - 4*Q9_9, 2*Q1_4 - 2*Q1_2 - 2*Q2_1 + 2*Q2_5 + 2*Q2_9 - 4*Q3_6 + 4*Q3_8 + 2*Q4_1 - 2*Q4_5 - 2*Q4_9 + 2*Q5_2 - 2*Q5_4 - 4*Q6_3 - 4*Q6_7 - 4*Q7_6 + 4*Q7_8 + 4*Q8_3 + 4*Q8_7 + 2*Q9_2 - 2*Q9_4, 4*Q2_6 - 2*Q1_7 - 2*Q1_3 + 4*Q2_8 - 2*Q3_1 + 2*Q3_5 - 2*Q3_9 + 4*Q4_6 + 4*Q4_8 + 2*Q5_3 + 2*Q5_7 + 4*Q6_2 + 4*Q6_4 - 2*Q7_1 + 2*Q7_5 - 2*Q7_9 + 4*Q8_2 + 4*Q8_4 - 2*Q9_3 - 2*Q9_7, 4*Q2_9 - 4*Q1_4 - 4*Q2_1 - 4*Q2_5 - 4*Q1_2 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 - 4*Q4_5 + 4*Q4_9 - 4*Q5_2 - 4*Q5_4 + 8*Q6_3 + 8*Q6_7 + 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 + 4*Q9_4, 8*Q3_3 - 8*Q2_2 + 8*Q4_4 - 8*Q6_6 - 8*Q7_7 + 8*Q8_8, 6*Q3_9 - 6*Q1_7 - 6*Q3_1 - 6*Q3_5 - 6*Q1_3 - 6*Q5_3 - 6*Q5_7 - 6*Q7_1 - 6*Q7_5 + 6*Q7_9 + 6*Q9_3 + 6*Q9_7, 4*Q1_6 - 4*Q1_8 - 8*Q2_3 - 8*Q2_7 - 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 + 4*Q5_6 - 4*Q5_8 + 4*Q6_1 + 4*Q6_5 - 4*Q6_9 - 8*Q7_2 + 8*Q7_4 - 4*Q8_1 - 4*Q8_5 + 4*Q8_9 - 4*Q9_6 + 4*Q9_8, 2*Q1_3 + 2*Q1_7 + 4*Q2_6 - 4*Q2_8 + 2*Q3_1 + 2*Q3_5 + 2*Q3_9 - 4*Q4_6 + 4*Q4_8 + 2*Q5_3 + 2*Q5_7 + 4*Q6_2 - 4*Q6_4 + 2*Q7_1 + 2*Q7_5 + 2*Q7_9 - 4*Q8_2 + 4*Q8_4 + 2*Q9_3 + 2*Q9_7, 2*Q5_6 - 2*Q1_8 - 2*Q1_6 + 2*Q5_8 - 2*Q6_1 + 2*Q6_5 - 2*Q6_9 - 2*Q8_1 + 2*Q8_5 - 2*Q8_9 - 2*Q9_6 - 2*Q9_8, 4*Q1_1 - 4*Q5_5 + 4*Q5_9 + 8*Q6_6 + 8*Q6_8 + 8*Q8_6 + 8*Q8_8 + 4*Q9_5 - 4*Q9_9, 2*Q1_2 - 2*Q1_4 + 2*Q2_1 - 2*Q2_5 + 2*Q2_9 + 4*Q3_6 + 4*Q3_8 - 2*Q4_1 + 2*Q4_5 - 2*Q4_9 - 2*Q5_2 + 2*Q5_4 + 4*Q6_3 - 4*Q6_7 - 4*Q7_6 - 4*Q7_8 + 4*Q8_3 - 4*Q8_7 + 2*Q9_2 - 2*Q9_4, 6*Q6_9 - 6*Q1_8 - 6*Q5_6 - 6*Q5_8 - 6*Q6_1 - 6*Q6_5 - 6*Q1_6 - 6*Q8_1 - 6*Q8_5 + 6*Q8_9 + 6*Q9_6 + 6*Q9_8, 4*Q1_7 - 4*Q1_3 - 8*Q2_6 - 8*Q2_8 - 4*Q3_1 - 4*Q3_5 + 4*Q3_9 + 8*Q4_6 + 8*Q4_8 - 4*Q5_3 + 4*Q5_7 - 8*Q6_2 + 8*Q6_4 + 4*Q7_1 + 4*Q7_5 - 4*Q7_9 - 8*Q8_2 + 8*Q8_4 + 4*Q9_3 - 4*Q9_7, 2*Q1_6 + 2*Q1_8 - 4*Q2_3 + 4*Q2_7 - 4*Q3_2 + 4*Q3_4 + 4*Q4_3 - 4*Q4_7 + 2*Q5_6 + 2*Q5_8 + 2*Q6_1 + 2*Q6_5 + 2*Q6_9 + 4*Q7_2 - 4*Q7_4 + 2*Q8_1 + 2*Q8_5 + 2*Q8_9 + 2*Q9_6 + 2*Q9_8, 4*Q1_1 + 4*Q1_5 - 4*Q1_9 + 4*Q5_1 + 4*Q5_5 - 4*Q5_9 - 4*Q9_1 - 4*Q9_5 + 4*Q9_9, 6*Q1_2 - 6*Q1_4 + 6*Q2_1 + 6*Q2_5 - 6*Q2_9 - 6*Q4_1 - 6*Q4_5 + 6*Q4_9 + 6*Q5_2 - 6*Q5_4 - 6*Q9_2 + 6*Q9_4, 8*Q2_2 - 4*Q1_5 - 4*Q1_1 - 8*Q2_4 - 8*Q4_2 + 8*Q4_4 - 4*Q5_1 - 4*Q5_5 + 4*Q9_9, 2*Q1_4 - 2*Q1_2 - 2*Q2_1 - 2*Q2_5 - 2*Q2_9 + 2*Q4_1 + 2*Q4_5 + 2*Q4_9 - 2*Q5_2 + 2*Q5_4 - 2*Q9_2 + 2*Q9_4];


%coefficients for hessian
h11 = [ 12*Q1_1 - 12*Q1_5 - 12*Q1_9 - 12*Q5_1 + 12*Q5_5 + 12*Q5_9 - 12*Q9_1 + 12*Q9_5 + 12*Q9_9, 12*Q1_2 + 12*Q1_4 + 12*Q2_1 - 12*Q2_5 - 12*Q2_9 + 12*Q4_1 - 12*Q4_5 - 12*Q4_9 - 12*Q5_2 - 12*Q5_4 - 12*Q9_2 - 12*Q9_4, 12*Q1_3 + 12*Q1_7 + 12*Q3_1 - 12*Q3_5 - 12*Q3_9 - 12*Q5_3 - 12*Q5_7 + 12*Q7_1 - 12*Q7_5 - 12*Q7_9 - 12*Q9_3 - 12*Q9_7, 12*Q1_8 - 12*Q1_6 + 12*Q5_6 - 12*Q5_8 - 12*Q6_1 + 12*Q6_5 + 12*Q6_9 + 12*Q8_1 - 12*Q8_5 - 12*Q8_9 + 12*Q9_6 - 12*Q9_8, 4*Q1_5 - 4*Q1_1 + 8*Q2_2 + 8*Q2_4 + 8*Q4_2 + 8*Q4_4 + 4*Q5_1 - 4*Q5_5 + 4*Q9_9, 4*Q1_6 + 4*Q1_8 + 8*Q2_3 + 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 - 4*Q5_6 - 4*Q5_8 + 4*Q6_1 - 4*Q6_5 - 4*Q6_9 + 8*Q7_2 + 8*Q7_4 + 4*Q8_1 - 4*Q8_5 - 4*Q8_9 - 4*Q9_6 - 4*Q9_8, 4*Q1_3 - 4*Q1_7 - 8*Q2_6 + 8*Q2_8 + 4*Q3_1 - 4*Q3_5 - 4*Q3_9 - 8*Q4_6 + 8*Q4_8 - 4*Q5_3 + 4*Q5_7 - 8*Q6_2 - 8*Q6_4 - 4*Q7_1 + 4*Q7_5 + 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 + 4*Q9_7, 4*Q1_9 - 4*Q1_1 + 8*Q3_3 + 8*Q3_7 + 4*Q5_5 + 8*Q7_3 + 8*Q7_7 + 4*Q9_1 - 4*Q9_9, 4*Q1_4 - 4*Q1_2 - 4*Q2_1 + 4*Q2_5 + 4*Q2_9 - 8*Q3_6 + 8*Q3_8 + 4*Q4_1 - 4*Q4_5 - 4*Q4_9 + 4*Q5_2 - 4*Q5_4 - 8*Q6_3 - 8*Q6_7 - 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 - 4*Q9_4, 4*Q1_1 - 4*Q5_5 - 4*Q5_9 + 8*Q6_6 - 8*Q6_8 - 8*Q8_6 + 8*Q8_8 - 4*Q9_5 - 4*Q9_9]; 
h12 = [ 6*Q1_2 + 6*Q1_4 + 6*Q2_1 - 6*Q2_5 - 6*Q2_9 + 6*Q4_1 - 6*Q4_5 - 6*Q4_9 - 6*Q5_2 - 6*Q5_4 - 6*Q9_2 - 6*Q9_4, 8*Q1_5 - 8*Q1_1 + 16*Q2_2 + 16*Q2_4 + 16*Q4_2 + 16*Q4_4 + 8*Q5_1 - 8*Q5_5 + 8*Q9_9, 4*Q1_6 + 4*Q1_8 + 8*Q2_3 + 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 - 4*Q5_6 - 4*Q5_8 + 4*Q6_1 - 4*Q6_5 - 4*Q6_9 + 8*Q7_2 + 8*Q7_4 + 4*Q8_1 - 4*Q8_5 - 4*Q8_9 - 4*Q9_6 - 4*Q9_8, 4*Q1_3 - 4*Q1_7 - 8*Q2_6 + 8*Q2_8 + 4*Q3_1 - 4*Q3_5 - 4*Q3_9 - 8*Q4_6 + 8*Q4_8 - 4*Q5_3 + 4*Q5_7 - 8*Q6_2 - 8*Q6_4 - 4*Q7_1 + 4*Q7_5 + 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 + 4*Q9_7, 6*Q2_5 - 6*Q1_4 - 6*Q2_1 - 6*Q1_2 - 6*Q2_9 - 6*Q4_1 + 6*Q4_5 - 6*Q4_9 + 6*Q5_2 + 6*Q5_4 - 6*Q9_2 - 6*Q9_4, 8*Q2_6 - 4*Q1_7 - 4*Q1_3 + 8*Q2_8 - 4*Q3_1 + 4*Q3_5 - 4*Q3_9 + 8*Q4_6 + 8*Q4_8 + 4*Q5_3 + 4*Q5_7 + 8*Q6_2 + 8*Q6_4 - 4*Q7_1 + 4*Q7_5 - 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 - 4*Q9_7, 4*Q1_6 - 4*Q1_8 + 8*Q2_3 - 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 - 8*Q4_7 - 4*Q5_6 + 4*Q5_8 + 4*Q6_1 - 4*Q6_5 + 4*Q6_9 - 8*Q7_2 - 8*Q7_4 - 4*Q8_1 + 4*Q8_5 - 4*Q8_9 + 4*Q9_6 - 4*Q9_8, 2*Q2_9 - 2*Q1_4 - 2*Q2_1 - 2*Q2_5 - 2*Q1_2 + 4*Q3_6 + 4*Q3_8 - 2*Q4_1 - 2*Q4_5 + 2*Q4_9 - 2*Q5_2 - 2*Q5_4 + 4*Q6_3 + 4*Q6_7 + 4*Q7_6 + 4*Q7_8 + 4*Q8_3 + 4*Q8_7 + 2*Q9_2 + 2*Q9_4, 8*Q3_3 - 8*Q2_2 + 8*Q4_4 - 8*Q6_6 - 8*Q7_7 + 8*Q8_8, 2*Q1_2 + 2*Q1_4 + 2*Q2_1 + 2*Q2_5 + 2*Q2_9 - 4*Q3_6 + 4*Q3_8 + 2*Q4_1 + 2*Q4_5 + 2*Q4_9 + 2*Q5_2 + 2*Q5_4 - 4*Q6_3 + 4*Q6_7 + 4*Q7_6 - 4*Q7_8 + 4*Q8_3 - 4*Q8_7 + 2*Q9_2 + 2*Q9_4]; 
h13 = [ 6*Q1_3 + 6*Q1_7 + 6*Q3_1 - 6*Q3_5 - 6*Q3_9 - 6*Q5_3 - 6*Q5_7 + 6*Q7_1 - 6*Q7_5 - 6*Q7_9 - 6*Q9_3 - 6*Q9_7, 4*Q1_6 + 4*Q1_8 + 8*Q2_3 + 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 - 4*Q5_6 - 4*Q5_8 + 4*Q6_1 - 4*Q6_5 - 4*Q6_9 + 8*Q7_2 + 8*Q7_4 + 4*Q8_1 - 4*Q8_5 - 4*Q8_9 - 4*Q9_6 - 4*Q9_8, 8*Q1_9 - 8*Q1_1 + 16*Q3_3 + 16*Q3_7 + 8*Q5_5 + 16*Q7_3 + 16*Q7_7 + 8*Q9_1 - 8*Q9_9, 4*Q1_4 - 4*Q1_2 - 4*Q2_1 + 4*Q2_5 + 4*Q2_9 - 8*Q3_6 + 8*Q3_8 + 4*Q4_1 - 4*Q4_5 - 4*Q4_9 + 4*Q5_2 - 4*Q5_4 - 8*Q6_3 - 8*Q6_7 - 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 - 4*Q9_4, 4*Q2_6 - 2*Q1_7 - 2*Q1_3 + 4*Q2_8 - 2*Q3_1 + 2*Q3_5 - 2*Q3_9 + 4*Q4_6 + 4*Q4_8 + 2*Q5_3 + 2*Q5_7 + 4*Q6_2 + 4*Q6_4 - 2*Q7_1 + 2*Q7_5 - 2*Q7_9 + 4*Q8_2 + 4*Q8_4 - 2*Q9_3 - 2*Q9_7, 4*Q2_9 - 4*Q1_4 - 4*Q2_1 - 4*Q2_5 - 4*Q1_2 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 - 4*Q4_5 + 4*Q4_9 - 4*Q5_2 - 4*Q5_4 + 8*Q6_3 + 8*Q6_7 + 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 + 4*Q9_4, 8*Q3_3 - 8*Q2_2 + 8*Q4_4 - 8*Q6_6 - 8*Q7_7 + 8*Q8_8, 6*Q3_9 - 6*Q1_7 - 6*Q3_1 - 6*Q3_5 - 6*Q1_3 - 6*Q5_3 - 6*Q5_7 - 6*Q7_1 - 6*Q7_5 + 6*Q7_9 + 6*Q9_3 + 6*Q9_7, 4*Q1_6 - 4*Q1_8 - 8*Q2_3 - 8*Q2_7 - 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 + 4*Q5_6 - 4*Q5_8 + 4*Q6_1 + 4*Q6_5 - 4*Q6_9 - 8*Q7_2 + 8*Q7_4 - 4*Q8_1 - 4*Q8_5 + 4*Q8_9 - 4*Q9_6 + 4*Q9_8, 2*Q1_3 + 2*Q1_7 + 4*Q2_6 - 4*Q2_8 + 2*Q3_1 + 2*Q3_5 + 2*Q3_9 - 4*Q4_6 + 4*Q4_8 + 2*Q5_3 + 2*Q5_7 + 4*Q6_2 - 4*Q6_4 + 2*Q7_1 + 2*Q7_5 + 2*Q7_9 - 4*Q8_2 + 4*Q8_4 + 2*Q9_3 + 2*Q9_7]; 

h22 = [ 4*Q1_5 - 4*Q1_1 + 8*Q2_2 + 8*Q2_4 + 8*Q4_2 + 8*Q4_4 + 4*Q5_1 - 4*Q5_5 + 4*Q9_9, 12*Q2_5 - 12*Q1_4 - 12*Q2_1 - 12*Q1_2 - 12*Q2_9 - 12*Q4_1 + 12*Q4_5 - 12*Q4_9 + 12*Q5_2 + 12*Q5_4 - 12*Q9_2 - 12*Q9_4, 8*Q2_6 - 4*Q1_7 - 4*Q1_3 + 8*Q2_8 - 4*Q3_1 + 4*Q3_5 - 4*Q3_9 + 8*Q4_6 + 8*Q4_8 + 4*Q5_3 + 4*Q5_7 + 8*Q6_2 + 8*Q6_4 - 4*Q7_1 + 4*Q7_5 - 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 - 4*Q9_7, 4*Q1_6 - 4*Q1_8 + 8*Q2_3 - 8*Q2_7 + 8*Q3_2 + 8*Q3_4 + 8*Q4_3 - 8*Q4_7 - 4*Q5_6 + 4*Q5_8 + 4*Q6_1 - 4*Q6_5 + 4*Q6_9 - 8*Q7_2 - 8*Q7_4 - 4*Q8_1 + 4*Q8_5 - 4*Q8_9 + 4*Q9_6 - 4*Q9_8, 12*Q1_1 - 12*Q1_5 + 12*Q1_9 - 12*Q5_1 + 12*Q5_5 - 12*Q5_9 + 12*Q9_1 - 12*Q9_5 + 12*Q9_9, 12*Q5_6 - 12*Q1_8 - 12*Q1_6 + 12*Q5_8 - 12*Q6_1 + 12*Q6_5 - 12*Q6_9 - 12*Q8_1 + 12*Q8_5 - 12*Q8_9 - 12*Q9_6 - 12*Q9_8, 12*Q1_7 - 12*Q1_3 - 12*Q3_1 + 12*Q3_5 - 12*Q3_9 + 12*Q5_3 - 12*Q5_7 + 12*Q7_1 - 12*Q7_5 + 12*Q7_9 - 12*Q9_3 + 12*Q9_7, 4*Q1_1 - 4*Q5_5 + 4*Q5_9 + 8*Q6_6 + 8*Q6_8 + 8*Q8_6 + 8*Q8_8 + 4*Q9_5 - 4*Q9_9, 4*Q1_2 - 4*Q1_4 + 4*Q2_1 - 4*Q2_5 + 4*Q2_9 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 + 4*Q4_5 - 4*Q4_9 - 4*Q5_2 + 4*Q5_4 + 8*Q6_3 - 8*Q6_7 - 8*Q7_6 - 8*Q7_8 + 8*Q8_3 - 8*Q8_7 + 4*Q9_2 - 4*Q9_4, 8*Q3_3 - 4*Q1_9 - 4*Q1_1 - 8*Q3_7 + 4*Q5_5 - 8*Q7_3 + 8*Q7_7 - 4*Q9_1 - 4*Q9_9];
h23 = [ 2*Q1_6 + 2*Q1_8 + 4*Q2_3 + 4*Q2_7 + 4*Q3_2 + 4*Q3_4 + 4*Q4_3 + 4*Q4_7 - 2*Q5_6 - 2*Q5_8 + 2*Q6_1 - 2*Q6_5 - 2*Q6_9 + 4*Q7_2 + 4*Q7_4 + 2*Q8_1 - 2*Q8_5 - 2*Q8_9 - 2*Q9_6 - 2*Q9_8, 8*Q2_6 - 4*Q1_7 - 4*Q1_3 + 8*Q2_8 - 4*Q3_1 + 4*Q3_5 - 4*Q3_9 + 8*Q4_6 + 8*Q4_8 + 4*Q5_3 + 4*Q5_7 + 8*Q6_2 + 8*Q6_4 - 4*Q7_1 + 4*Q7_5 - 4*Q7_9 + 8*Q8_2 + 8*Q8_4 - 4*Q9_3 - 4*Q9_7, 4*Q2_9 - 4*Q1_4 - 4*Q2_1 - 4*Q2_5 - 4*Q1_2 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 - 4*Q4_5 + 4*Q4_9 - 4*Q5_2 - 4*Q5_4 + 8*Q6_3 + 8*Q6_7 + 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 + 4*Q9_4, 8*Q3_3 - 8*Q2_2 + 8*Q4_4 - 8*Q6_6 - 8*Q7_7 + 8*Q8_8, 6*Q5_6 - 6*Q1_8 - 6*Q1_6 + 6*Q5_8 - 6*Q6_1 + 6*Q6_5 - 6*Q6_9 - 6*Q8_1 + 6*Q8_5 - 6*Q8_9 - 6*Q9_6 - 6*Q9_8, 8*Q1_1 - 8*Q5_5 + 8*Q5_9 + 16*Q6_6 + 16*Q6_8 + 16*Q8_6 + 16*Q8_8 + 8*Q9_5 - 8*Q9_9, 4*Q1_2 - 4*Q1_4 + 4*Q2_1 - 4*Q2_5 + 4*Q2_9 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 + 4*Q4_5 - 4*Q4_9 - 4*Q5_2 + 4*Q5_4 + 8*Q6_3 - 8*Q6_7 - 8*Q7_6 - 8*Q7_8 + 8*Q8_3 - 8*Q8_7 + 4*Q9_2 - 4*Q9_4, 6*Q6_9 - 6*Q1_8 - 6*Q5_6 - 6*Q5_8 - 6*Q6_1 - 6*Q6_5 - 6*Q1_6 - 6*Q8_1 - 6*Q8_5 + 6*Q8_9 + 6*Q9_6 + 6*Q9_8, 4*Q1_7 - 4*Q1_3 - 8*Q2_6 - 8*Q2_8 - 4*Q3_1 - 4*Q3_5 + 4*Q3_9 + 8*Q4_6 + 8*Q4_8 - 4*Q5_3 + 4*Q5_7 - 8*Q6_2 + 8*Q6_4 + 4*Q7_1 + 4*Q7_5 - 4*Q7_9 - 8*Q8_2 + 8*Q8_4 + 4*Q9_3 - 4*Q9_7, 2*Q1_6 + 2*Q1_8 - 4*Q2_3 + 4*Q2_7 - 4*Q3_2 + 4*Q3_4 + 4*Q4_3 - 4*Q4_7 + 2*Q5_6 + 2*Q5_8 + 2*Q6_1 + 2*Q6_5 + 2*Q6_9 + 4*Q7_2 - 4*Q7_4 + 2*Q8_1 + 2*Q8_5 + 2*Q8_9 + 2*Q9_6 + 2*Q9_8];

h33 = [ 4*Q1_9 - 4*Q1_1 + 8*Q3_3 + 8*Q3_7 + 4*Q5_5 + 8*Q7_3 + 8*Q7_7 + 4*Q9_1 - 4*Q9_9, 4*Q2_9 - 4*Q1_4 - 4*Q2_1 - 4*Q2_5 - 4*Q1_2 + 8*Q3_6 + 8*Q3_8 - 4*Q4_1 - 4*Q4_5 + 4*Q4_9 - 4*Q5_2 - 4*Q5_4 + 8*Q6_3 + 8*Q6_7 + 8*Q7_6 + 8*Q7_8 + 8*Q8_3 + 8*Q8_7 + 4*Q9_2 + 4*Q9_4, 12*Q3_9 - 12*Q1_7 - 12*Q3_1 - 12*Q3_5 - 12*Q1_3 - 12*Q5_3 - 12*Q5_7 - 12*Q7_1 - 12*Q7_5 + 12*Q7_9 + 12*Q9_3 + 12*Q9_7, 4*Q1_6 - 4*Q1_8 - 8*Q2_3 - 8*Q2_7 - 8*Q3_2 + 8*Q3_4 + 8*Q4_3 + 8*Q4_7 + 4*Q5_6 - 4*Q5_8 + 4*Q6_1 + 4*Q6_5 - 4*Q6_9 - 8*Q7_2 + 8*Q7_4 - 4*Q8_1 - 4*Q8_5 + 4*Q8_9 - 4*Q9_6 + 4*Q9_8, 4*Q1_1 - 4*Q5_5 + 4*Q5_9 + 8*Q6_6 + 8*Q6_8 + 8*Q8_6 + 8*Q8_8 + 4*Q9_5 - 4*Q9_9, 12*Q6_9 - 12*Q1_8 - 12*Q5_6 - 12*Q5_8 - 12*Q6_1 - 12*Q6_5 - 12*Q1_6 - 12*Q8_1 - 12*Q8_5 + 12*Q8_9 + 12*Q9_6 + 12*Q9_8, 4*Q1_7 - 4*Q1_3 - 8*Q2_6 - 8*Q2_8 - 4*Q3_1 - 4*Q3_5 + 4*Q3_9 + 8*Q4_6 + 8*Q4_8 - 4*Q5_3 + 4*Q5_7 - 8*Q6_2 + 8*Q6_4 + 4*Q7_1 + 4*Q7_5 - 4*Q7_9 - 8*Q8_2 + 8*Q8_4 + 4*Q9_3 - 4*Q9_7, 12*Q1_1 + 12*Q1_5 - 12*Q1_9 + 12*Q5_1 + 12*Q5_5 - 12*Q5_9 - 12*Q9_1 - 12*Q9_5 + 12*Q9_9, 12*Q1_2 - 12*Q1_4 + 12*Q2_1 + 12*Q2_5 - 12*Q2_9 - 12*Q4_1 - 12*Q4_5 + 12*Q4_9 + 12*Q5_2 - 12*Q5_4 - 12*Q9_2 + 12*Q9_4, 8*Q2_2 - 4*Q1_5 - 4*Q1_1 - 8*Q2_4 - 8*Q4_2 + 8*Q4_4 - 4*Q5_1 - 4*Q5_5 + 4*Q9_9]; 


% strat the iteration
% s0t = permute(s0,[2 1 3]);
% Q = repmat(Q,1,1,size(s0,3));
s = s0;
cost_pre = cost_function_val(s,K);
eye = [1 0 0; 0 1 0; 0 0 1];
for i = 1 : iter 
    s1 = s(1);  s2 = s(2);  s3 = s(3);
    vec = [ s1^2, s1*s2, s1*s3, s1, s2^2, s2*s3, s2, s3^2, s3, 1]';
    vec3 = [s1^3, s1^2*s2, s1^2*s3, s1^2, s1*s2^2, s1*s2*s3, s1*s2, s1*s3^2, s1*s3, s1, s2^3, s2^2*s3, s2^2, s2*s3^2, s2*s3, s2, s3^3, s3^2, s3, 1]';
    %gradient
    G = [coeff1*vec3; coeff2*vec3; coeff3*vec3];
    
    %hessian matrix
    H = [h11*vec h12*vec h13*vec; 
     h12*vec h22*vec h23*vec; 
     h13*vec h23*vec h33*vec ];
     if min(eig(H)) < -1e-3
        s = zeros(3,1);
        return;
     end
 
    d = ( H + lambda * eye )\G;
    s = s - d;
    cost = cost_function_val(s,K);
    
    if cost < cost_pre
        lambda = lambda/10;
        if lambda < 1
            lambda = 1;
        end
    else
        lambda = lambda*10;
    end
    cost_pre = cost;          
end

end

