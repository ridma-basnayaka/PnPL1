function [N, D] = buildNDMatrices(points3D, points2D, lines3D, lines2D)
% Build N and D matrices from point and line constraints

N = [];
D = [];

% Process point correspondences
for i = 1:size(points3D, 2)
    p = points3D(:, i);
    b = points2D(1:2, i); % Extract u,v (ignore homogeneous 1)
    u = b(1); v = b(2);
    
    % First equation: -p' for r1, 0 for r2, u*p' for r3
    D_i1 = [-p', zeros(1,3), u*p'];
    N_i1 = [-1, 0, u];
    
    % Second equation: 0 for r1, -p' for r2, v*p' for r3
    D_i2 = [zeros(1,3), -p', v*p'];
    N_i2 = [0, -1, v];
    
    D = [D; D_i1; D_i2];
    N = [N; N_i1; N_i2];
end

% Process line correspondences
for i = 1:size(lines3D.directions, 2)
    d = lines3D.directions(:, i);
    p = lines3D.points(:, i);
    n = lines2D.normals(:, i);
    
    % Line orientation constraint (1 equation)
    D_orient = kron(d', n');
    N_orient = zeros(1, 3);
    
    % Line position constraint (1 equation)
    D_pos = kron(p', n');
    N_pos = n';
    
    D = [D; D_orient; D_pos];
    N = [N; N_orient; N_pos];
end
end