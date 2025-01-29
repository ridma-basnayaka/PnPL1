function R = quaternionRotation(s)
% Compute rotation matrix from quaternion parameters s = [s1, s2, s3, s4]

s1 = s(1); s2 = s(2); s3 = s(3); s4 = s(4);
Q = s1^2 + s2^2 + s3^2 + s4^2;

R = [ (s1^2 + s2^2 - s3^2 - s4^2), 2*(s2*s3 - s1*s4), 2*(s2*s4 + s1*s3);
      2*(s2*s3 + s1*s4), (s1^2 - s2^2 + s3^2 - s4^2), 2*(s3*s4 - s1*s2);
      2*(s2*s4 - s1*s3), 2*(s3*s4 + s1*s2), (s1^2 - s2^2 - s3^2 + s4^2) ] / Q;
end