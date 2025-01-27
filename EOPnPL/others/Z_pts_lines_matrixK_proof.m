Wp = sym('x',[3,8],'real');
Wl = sym('y',[3,8],'real');
W = [Wp;Wl];
f1 = W'*W;
f2 = Wp'*Wp + Wl'*Wl;
disp(f1 - f2);

Vp = sym('m',[3,8],'real');
Vl = sym('n',[3,8],'real');
V = [Vp;Vl];
f3 = W'*V; 
f4 = Wp'*Vp + Wl'*Vl;
disp(f3 - f4);

