f = xrand(1,1,[200 2000]);
width = 640;
height = 480;
A = [f, 0, width/2; 0 0 height/2; 0 0 1];

% generate 3d coordinates in camera space
Xc= [xrand(1,10,[-2 2]); xrand(1,10,[-2 2]); xrand(1,10,[4 8])];
t= mean(Xc,2);
R= rodrigues(randn(3,1));
P = A * [R,t];

[Ad,Rd,td] = P2Rt(P);