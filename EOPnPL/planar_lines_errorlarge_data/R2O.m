function [v] = R2O(R)

    x = atan2(R(3,2),R(3,3));
    y = atan2(-R(3,1),sqrt(R(3,2)^2+R(3,3)^2));
    z = atan2(R(2,1),R(1,1));
    degree_x = x / (2*pi) * 180;
    degree_y = y / (2*pi) * 180;
    degree_z = z / (2*pi) * 180;
    v = [x,y,z,degree_x,degree_y,degree_z];
end