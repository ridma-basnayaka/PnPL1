function [] = savetxt(fname,data)
%将数据存储为txt文本
%一个数据换一行
address = 'D:\实现\EOPnP\EOPnP\';
format = '.txt';
addr_name_for = [address,fname,format];
fid = fopen(addr_name_for,'w+');
[m,n] = size(data);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%g\n',data(i,j));
    end
end
fclose(fid);

