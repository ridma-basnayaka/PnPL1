function Johnny_plot_total_comparison
load('Johnny_total_comparison.mat');

close all;
addpath('util\');
h = sp_position2();
set(gcf,'Units','normal');
set(gcf, 'PaperPositionMode', 'auto');

% subplot(1,5,1);
sp_format2(2);
XLabel = 'Number of constraints';
Xarg = ns;
ws = ones(length(method_list), 1);
yrange = [80 200];
[mnames, p] = xdraw_main(Xarg,yrange,method_list,'accur_med_r','Median Rotation Ratio',XLabel,'Ratio (%)', ws);
correct_margin();
% {'Median Rotation Error'; 'Ratio to OPnPL'}
% subplot(1,5,2);
sp_format2(1);
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'accur_mean_r','Mean Rotation Ratio',...
    XLabel,'Ratio (%)',ws);
correct_margin();

% subplot(1,5,3);
sp_format2(4);
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'accur_med_t','Median Translation Ratio',...
    XLabel,'Ratio (%)',ws);
correct_margin();

% subplot(1, 5, 4);
sp_format2(3);
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'accur_mean_t','Mean Translation Ratio',...
    XLabel,'Ratio (%)',ws);
correct_margin();

% subplot(1, 5, 5);
sp_format2(5);
yrange = [0 120];
ws = ones(length(method_list), 1);
xdraw_main(Xarg,yrange,method_list,'efficiency_ratio','Mean Runtime Ratio',...
    XLabel,'Ratio (%)',ws);
correct_margin();

% hL = legend(mnames, 'Orientation','horizontal', 'Position', [0.2 0.9 0.6 0.1]);
hL = legend(mnames, 'Orientation','vertical', 'Position', [0.88 0.2 0.1 0.6]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
end