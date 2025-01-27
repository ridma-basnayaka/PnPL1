function [R, t] = AlgLS_lx(xs,xe,Ps,Pe)
%This function follows the AlgLS algorithm proposed by Faraz and Stergios: 
% Globally optimal pose estimation from line corrsepondences.
%李鑫修改了接口和输出
%输入像点起点和终点，物点起点和终点
%输出 从世界系到像机系的变换，R和t 与LuCP的一致，与我常用的pnp程序一致
n=length(xs);
xs1=[xs;ones(1,n)];
xe1=[xe;ones(1,n)];
Pw=Ps;
Vw=Pe-Ps;
for i=1:n
    Vw(:,i)=Vw(:,i)/norm(Vw(:,i));
end
[rot_cw, pos_cw] = AlgLS(xs1,xe1,Vw,Pw);
R=rot_cw';
t=-R*pos_cw;
end
