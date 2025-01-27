function [R, t] = AlgLS_lx(xs,xe,Ps,Pe)
%This function follows the AlgLS algorithm proposed by Faraz and Stergios: 
% Globally optimal pose estimation from line corrsepondences.
%�����޸��˽ӿں����
%������������յ㣬��������յ�
%��� ������ϵ�����ϵ�ı任��R��t ��LuCP��һ�£����ҳ��õ�pnp����һ��
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
