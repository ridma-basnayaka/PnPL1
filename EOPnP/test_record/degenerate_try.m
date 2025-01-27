A =     [   1000        0           320;
            0           1000        240;
            0           0           1   ];
number = 10;
Xw = zeros(3,number);
for i = 1:number
    x = -2 + 4*rand();
    y = -2 + 4*rand();
    z = 4  + 4*rand();
    Xw(:,i) = [x;y;z];
end
t1 = -5 + 10*rand();
t2 = -5 + 10*rand();
t3 = 0 + 15*rand(); %t3不能小于零，即不能在摄像机后面
t = [t1;t2;t3];

%退化情况
Rd = CGR_degenerate();

Ximg_hd = A * ( Rd * Xw + t );
Ximgd(1,:) = Ximg_hd(1,:) ./ Ximg_hd(3,:);
Ximgd(2,:) = Ximg_hd(2,:) ./ Ximg_hd(3,:);
Ximgd = vpa(Ximgd);
[Red,ted,errfd,errd,s0,s] = EOPnP(Ximgd,Xw,A);