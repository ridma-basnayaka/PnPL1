function [points3D, points2D, lines3D, lines2D] = generateSyntheticData()
% Generate synthetic camera pose estimation problem

% Ground truth pose
theta = pi/6;
R_gt = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
t_gt = [0.5; -0.2; 1.0];

% Generate 3D points
n_p = 6;
points3D = 5*randn(3, n_p);

% Project points (convert to normalized coordinates)
points2D = R_gt * points3D + t_gt;
points2D = points2D ./ points2D(3,:); % Normalize

% Generate 3D lines
n_l = 4;
lines3D.directions = randn(3, n_l);
lines3D.points = 5*randn(3, n_l);

% Project lines
lines2D.normals = zeros(3, n_l);
for i = 1:n_l
    d = lines3D.directions(:,i);
    p = lines3D.points(:,i);
    
    % Projected line normal (correct formula)
    line_dir_cam = R_gt * d;
    line_pt_cam = R_gt * p + t_gt;
    n = cross(line_dir_cam, line_pt_cam);
    lines2D.normals(:,i) = n/norm(n);
end
end