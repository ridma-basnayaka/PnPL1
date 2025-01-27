addpath others
R_temp1 = zeros(3,3,0); t_temp1 = zeros(3,0);
Rn = cat(3, rotx(pi/2), roty(pi/2), rotz(pi/2));
tn = mean([XXw Xs Xe],2);

cost = inf;


for i = 1:3
    % Make a random rotation
    XXw_n = Rn(:,:,i) * (XXw - repmat(tn, 1, size(XXw,2)));
    Xs_n = Rn(:,:,i) * (Xs - repmat(tn, 1, size(Xs,2)));
    Xe_n = Rn(:,:,i) * (Xe - repmat(tn, 1, size(Xe,2)));
    
    [Ri, ti] = EOPnPL(XXw_n,xxn,xs,xe,Xs_n,Xe_n);
     
    R_temp1 = cat(3,R_temp1, Ri);
    t_temp1 = [t_temp1 ti];
end

[minInd,~,~,~] = findBestRTReproj(R_temp1, t_temp1, xs, xe, Xs, Xe, XXw, xxn);
R_rob = R_temp1(:,:,minInd);
t_rob = t_temp1(:,minInd);

tempy_rob = cal_pose_err([R_rob t_rob],[R t]);
y_rob = tempy_rob;
disp(y_rob);


function r = rotx(t)
ct = cos(t);
st = sin(t);
r =    [1	0	0;
    0	ct	-st;
    0	st	ct];
end

function r = roty(t)
% roty: rotation about y-axi-
ct = cos(t);
st = sin(t);
r =    [ct	0	st;
    0	1	0;
    -st	0	ct];
end

function r = rotz(t)
% rotz: rotation about z-axis

ct = cos(t);
st = sin(t);
r = [ct	-st	0
    st	ct	0
    0	0	1];

end
