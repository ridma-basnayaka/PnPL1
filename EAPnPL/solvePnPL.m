function [R_best, t_best] = solvePnPL(N, D)
% Solve for camera pose using PnPL constraints

% Construct matrix A (Equation 14)
I = eye(size(N,1));
A = (N * ((N' * N) \ N') - I) * D;

% Build coefficient matrix C
C = A' * A;

% Eigen decomposition for quaternion parameters
[V, ~] = eig(C);
quats = V(:,1:4); % First 4 eigenvectors

% Evaluate all possible solutions
min_error = inf;
for k = 1:size(quats,2)
    s = quats(:,k);
    s = s/norm(s);
    
    % Convert to rotation matrix
    R = quaternionRotation(s(1:4));
    
    % Calculate translation
    t = (N'*N)\(N'*D*R(:));
    
    % Calculate residual error
    error = norm(A*R(:))^2;
    
    if error < min_error
        min_error = error;
        R_best = R;
        t_best = t;
    end
end
end