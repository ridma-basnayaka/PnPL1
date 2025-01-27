function plot_timings
    load('lines_pts_number_timing');
    Xarg = 2*ns;
%     Xarg = ns;
    for i = 1:length(method_list)
%         method_list(i).avg_t_acc1 = method_list(i).avg_t_1 + method_list(i).avg_t_2;
%         method_list(i).avg_t_acc2 = method_list(i).avg_t_2 + method_list(i).avg_t_4;
        method_list(i).avg_t_acc1 = method_list(i).avg_t_1*1000;
    end
    

    w= 300;
    h= 300;
    XLabel = 'Number of constraints';
    close all;    
    figure('color','w','position',[0,0,w,h]);
    yrange = [0 35];    
    xdrawgraph(Xarg,yrange,method_list,'avg_t_acc1','Mean processing time',...
    XLabel,'Computational Time(milliseconds)', ones(size(method_list,2),1));

%     figure('color','w','position',[0,0,w,h]);
%     yrange = [0 4e-2];    
%     xdrawgraph(Xarg,yrange,method_list,'avg_t_acc2','Mean solving time',...
%     XLabel,'Time, s.', ones(2,1));


end