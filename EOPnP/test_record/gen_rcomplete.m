r = sym('r',[6,2],'real');
%恢复完整的r
c1 = r(1:2:5,:);    %第一列
c2 = r(2:2:6,:);    %第二列
n = size(r,2);
rcomplete = sym('rcomplete',[9,n],'real');
for i = 1 : n
    c3 = xcross(c1(:,i),c2(:,i));
    rcomplete(1:2,i) = r(1:2,i);
    rcomplete(4:5,i) = r(3:4,i);
    rcomplete(7:8,i) = r(5:6,i);
    rcomplete(3:3:9,i) = c3;
end
disp(r);
disp(rcomplete);

disp([  r(3,:).*r(6,:) - r(4,:).*r(5,:);
        r(2,:).*r(5,:) - r(1,:).*r(6,:);
        r(1,:).*r(4,:) - r(2,:).*r(3,:)]);

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
end