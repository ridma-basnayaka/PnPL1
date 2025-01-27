n = 4;
u = sym('u',[1,n],'real');
v = sym('v',[1,n],'real');
x = sym('x',[1,n],'real');
y = sym('y',[1,n],'real');
z = sym('z',[1,n],'real');
Ximgn = [u;v];
Xw = [x;y;z];
w1 = repmat(eye(2),n,1);
w2 = (-1) *reshape(Ximgn,2*n,1);
w = [w1 w2];
wtw = w'*w;
disp('WtW = ')
disp(wtw);  %观察易得规律
wtw_inv = inv(wtw);
wtw_inv = simplify(wtw_inv);    %3*3矩阵

%固定 wtw_inv，继续分解
m = sym('m',[3,3],'real');
m1 = m(1:2,1:2);
m2 = m(1:2,3);
m3 = m(3,1:2);
m4 = m(3,3);    %m = [m1 m2; m3 m4]
wm = simplify( w * m );
wm_decompose = [w1*m1+w2*m3, w1*m2+w2*m4];
simplify(wm - wm_decompose );

%分解W*(WT*W)^-1*WT
WWT = w * m * w';
wwt = w1 * m1 * w1' + w2 * m3 * w1' + w1 * m2 * w2' + w2 * m4 * w2';
simplify( WWT - wwt );



wmwt = sym('wmwt',[2*n,2*n],'real');
for i = 0:2*n-1
   for j = 0:2*n-1
      wmwt(i+1,j+1) =  m1(mod(i,2)+1,mod(j,2)+1) + (-1)*Ximgn(mod(i,2)+1,floor(i/2)+1)*m3(1,mod(j,2)+1) + (-1)*m2(mod(i,2)+1,1)*Ximgn(mod(j,2)+1,floor(j/2)+1) +...
          m4*Ximgn(mod(i,2)+1,floor(i/2)+1)*Ximgn(mod(j,2)+1,floor(j/2)+1);
   end
end

 simplify( WWT - wmwt )


for i = 0:n-1
    V(2*i+1,1:3) = -Xw(:,i+1)';
    V(2*i+1,4:6) = 0;
    V(2*i+2,1:3) = 0;
    V(2*i+2,4:6) = -Xw(:,i+1)';
    V(2*i+1,7:9) = u(i+1) * Xw(:,i+1)';
    V(2*i+2,7:9) = v(i+1) * Xw(:,i+1)';
end

W_new = sym('W_new',[2*n,2*n],'real');
X = -Xw(1,:);
Y = -Xw(2,:);
Z = -Xw(3,:);
XU = -X.*Ximgn(1,:);
XV = -X.*Ximgn(2,:);
YU = -Y.*Ximgn(1,:);
YV = -Y.*Ximgn(2,:);
ZU = -Z.*Ximgn(1,:);
ZV = -Z.*Ximgn(2,:);
W_new1 = W_new(:,1:2:2*n-1);
W_new2 = W_new(:,2:2:2*n);
f1 = W_new*V - V;
for i = 1:2*n
    for j = 1:9
        f2(i,1) = dot(W_new1(i,:),X) - V(i,1);
        f2(i,2) = dot(W_new1(i,:),Y) - V(i,2);
        f2(i,3) = dot(W_new1(i,:),Z) - V(i,3);
        f2(i,4) = dot(W_new2(i,:),X) - V(i,4);
        f2(i,5) = dot(W_new2(i,:),Y) - V(i,5);
        f2(i,6) = dot(W_new2(i,:),Z) - V(i,6);
        f2(i,7) = dot(W_new1(i,:),XU) + dot(W_new2(i,:),XV) - V(i,7);
        f2(i,8) = dot(W_new1(i,:),YU) + dot(W_new2(i,:),YV) - V(i,8);
        f2(i,9) = dot(W_new1(i,:),ZU) + dot(W_new2(i,:),ZV) - V(i,9);
    end
end
f1-f2



