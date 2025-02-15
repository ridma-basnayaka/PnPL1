function [r] = Planar_case_1_z(v1)

v11 = v1(1);
v12 = v1(2);
v13 = v1(3);
v14 = v1(4);
v15 = v1(5);
v16 = v1(6);

coeffs_a1 = zeros(1,7);
coeffs_a1(7:-2:1) = [ 4*(v11*v14 - v12*v13)^3*(2*v11*v14 - 2*v12*v13) + 4*(v11*v16 - v12*v15)^3*(2*v11*v16 - 2*v12*v15) + 4*(v13*v16 - v14*v15)^3*(2*v13*v16 - 2*v14*v15) + 2*(v11*v14 - v12*v13)*(v11*v16 - v12*v15)*((v11*v14 - v12*v13)*(2*v11*v16 - 2*v12*v15) + (2*v11*v14 - 2*v12*v13)*(v11*v16 - v12*v15)) + 2*(v11*v14 - v12*v13)*(v13*v16 - v14*v15)*((v11*v14 - v12*v13)*(2*v13*v16 - 2*v14*v15) + (2*v11*v14 - 2*v12*v13)*(v13*v16 - v14*v15)) + 2*(v11*v16 - v12*v15)*(v13*v16 - v14*v15)*((v11*v16 - v12*v15)*(2*v13*v16 - 2*v14*v15) + (2*v11*v16 - 2*v12*v15)*(v13*v16 - v14*v15)), 2*(v11*v15 + v12*v16)*((v11*v14 - v12*v13)*(2*v13*v16 - 2*v14*v15) + (2*v11*v14 - 2*v12*v13)*(v13*v16 - v14*v15)) - 2*(v13*v15 + v14*v16)*((v11*v14 - v12*v13)*(2*v11*v16 - 2*v12*v15) + (2*v11*v14 - 2*v12*v13)*(v11*v16 - v12*v15)) - 2*(v11*v13 + v12*v14)*((v11*v16 - v12*v15)*(2*v13*v16 - 2*v14*v15) + (2*v11*v16 - 2*v12*v15)*(v13*v16 - v14*v15)) + 2*(2*v11^2 + 2*v12^2)*(v13*v16 - v14*v15)^2 + 2*(2*v13^2 + 2*v14^2)*(v11*v16 - v12*v15)^2 + 2*(2*v15^2 + 2*v16^2)*(v11*v14 - v12*v13)^2 - 2*(v11*v14 - v12*v13)*(v11*v16 - v12*v15)*(2*v13*v15 + 2*v14*v16) + 2*(v11*v14 - v12*v13)*(2*v11*v15 + 2*v12*v16)*(v13*v16 - v14*v15) - 2*(2*v11*v13 + 2*v12*v14)*(v11*v16 - v12*v15)*(v13*v16 - v14*v15) + 4*(v15^2 + v16^2)*(v11*v14 - v12*v13)*(2*v11*v14 - 2*v12*v13) + 4*(v13^2 + v14^2)*(v11*v16 - v12*v15)*(2*v11*v16 - 2*v12*v15) + 4*(v11^2 + v12^2)*(v13*v16 - v14*v15)*(2*v13*v16 - 2*v14*v15), 2*(v11^2 + v12^2)*(2*v11^2 + 2*v12^2) + 2*(v13^2 + v14^2)*(2*v13^2 + 2*v14^2) + 2*(v15^2 + v16^2)*(2*v15^2 + 2*v16^2) + 2*(v11*v13 + v12*v14)*(2*v11*v13 + 2*v12*v14) - 4*(v11*v14 - v12*v13)*(2*v11*v14 - 2*v12*v13) + 2*(v11*v15 + v12*v16)*(2*v11*v15 + 2*v12*v16) - 4*(v11*v16 - v12*v15)*(2*v11*v16 - 2*v12*v15) + 2*(v13*v15 + v14*v16)*(2*v13*v15 + 2*v14*v16) - 4*(v13*v16 - v14*v15)*(2*v13*v16 - 2*v14*v15), - 4*v11^2 - 4*v12^2 - 4*v13^2 - 4*v14^2 - 4*v15^2 - 4*v16^2];
solution_a1 = roots(coeffs_a1);
a1 = solution_a1(abs(imag(solution_a1)) < 1e-10);    % 取实数
r = v1 * a1';

end

