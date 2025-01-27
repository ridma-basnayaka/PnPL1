% parameter setting 1(translation = 0) equations check
% 
% 
% (R_true(2,:) * Xw + t_true12(2)) - Ximg_h(2,:).*(R_true(3,:) * Xw + t_true12(3))
% (R_true(1,:) * Xw + t_true12(1)) - Ximg_h(1,:).*(R_true(3,:) * Xw + t_true12(3))
% Ximg_h = A\[Ximg;ones(1,size(Ximg,2))]
% (R_true(2,:) * Xw + t_true12(2)) ./ (R_true(3,:) * Xw + t_true12(3))
% (R_true(1,:) * Xw + t_true12(1)) ./ (R_true(3,:) * Xw + t_true12(3))
% (R_bar(2,:) * Xw + (1+n'*n)*t_true12(2)) - Ximg_h(2,:).*(R_bar(3,:) * Xw + (1+n'*n)*t_true12(3))
% (R_bar(1,:) * Xw + (1+n'*n)*t_true12(1)) - Ximg_h(1,:).*(R_bar(3,:) * Xw + (1+n'*n)*t_true12(3))
% T12 = (1+n'*n).*t_true12
% (R_bar(2,:) * Xw + T12(2)) - Ximg_h(2,:).*(R_bar(3,:) * Xw + T12(3))
% (R_bar(1,:) * Xw + T12(1)) - Ximg_h(1,:).*(R_bar(3,:) * Xw + T12(3))
% 
% 
%  parameter setting 3(translation != 0) equations check
% Ximg_h3 = A\[Ximg3;ones(1,size(Ximg3,2))]
% (R_true(2,:) * Xw3 + t_true3(2)) ./ (R_true(3,:) * Xw3 + t_true3(3))
% (R_true(1,:) * Xw3 + t_true3(1)) ./ (R_true(3,:) * Xw3 + t_true3(3))
% T3 = (1+n'*n).*t_true3
% (R_bar(2,:) * Xw3 + T3(2)) - Ximg_h3(2,:).*(R_bar(3,:) * Xw3 + T3(3))
% (R_bar(1,:) * Xw3 + T3(1)) - Ximg_h3(1,:).*(R_bar(3,:) * Xw3 + T3(3))
% 
% 
% 
% 
% 
% camera2 (translation = 0) check
% Ximg_h = A\[Ximg;ones(1,size(Ximg,2))]
% T = (1+s'*s).*t_true
% (R_true(1,:) * Xw + t_true(1)) ./ (R_true(3,:) * Xw + t_true(3))
% (R_true(2,:) * Xw + t_true(2)) ./ (R_true(3,:) * Xw + t_true(3))
% (R_bar(1,:) * Xw + T(1)) - Ximg_h(1,:).*(R_bar(3,:) * Xw + T(3))
% (R_bar(2,:) * Xw + T(2)) - Ximg_h(2,:).*(R_bar(3,:) * Xw + T(3))
% 
% 
% 检查罗德里格斯表达式
% 
% eye(3) + (1-cos(theta)).*k_x*k_x + sin(theta).*k_x
% cos(theta).*eye(3)+(1-cos(theta)).*k*k'+sin(theta)*k_x