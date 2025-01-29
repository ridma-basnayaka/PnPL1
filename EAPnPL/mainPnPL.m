% PnPL Camera Pose Estimation Main Script
clear; clc; close all;

% Generate synthetic data
[points3D, points2D, lines3D, lines2D] = generateSyntheticData();

% Build constraint matrices
[N, D] = buildNDMatrices(points3D, points2D(1:2,:), lines3D, lines2D);

% Solve PnPL problem
[R_est, t_est] = solvePnPL(N, D);

% Display results
disp('Estimated Rotation:');
disp(R_est);
disp('Estimated Translation:');
disp(t_est);

% Calculate angular error
R_gt = [cos(pi/6) 0 sin(pi/6); 0 1 0; -sin(pi/6) 0 cos(pi/6)];
angular_error = acosd((trace(R_gt'*R_est) - 1)/2);
fprintf('\nRotation Error: %.2f degrees\n', angular_error);

% Calculate translation error
translation_error = norm(t_est - [0.5; -0.2; 1.0]);
fprintf('Translation Error: %.4f m\n', translation_error); 
