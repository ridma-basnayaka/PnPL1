clear; 
clc;
close all;
warning off;

prepare_paths();

% experimental parameters

%focal length in pixels
fpix = 800;
%varying numbers of points and lines
ns = [10:100:1010];
%length of a line segment in 3D 
segLenBig = 3;
%scale of shift between line segments projected and present in the model
segLenShiftScale = 0.5*segLenBig;
%noise in pix, added to point projections, std. dev.
%nll = 1;
nll = 2;
%noise in pix, added to line segment endpoints' projections, std. dev.
% nlp = 1;
nlp = 2;
%number of random trials per iteration
num = 100;

%methods chosen for comparison
method_names = {'OPnPL', 'EOPnPL','EOPnPL-BJ-0','EPnPL'};
% method_names = {'OPnPL'};
is_planar = 0;
method_list = model.setup_methods(method_names, num, length(ns), is_planar);


%experiments
for i= 1:length(ns)
            
    fprintf('ns = %d: ',ns(i));

    index_fail = [];
    
    for k = 1:length(method_list)
        method_list(k).fails(i) = 0;
    end   
    
    nlines = ns(i);
    npt = ns(i);
    for j= 1:num
        
        % camera's parameters

        [R, t, XLTw, XXw, xxn, xs,xe,Xs,Xe] = model.setup_3d_scene(npt, nlines, segLenBig, segLenShiftScale, nll, nlp, fpix);

        for k = 1:length(method_list)
            
            % pose estimation
            [R1, t1, is_fail, s] = model.evaluate_method_timing(method_list(k), Xs, Xe, xs, xe, XXw, xxn, XLTw);      
            if (size(R1,1) < 3 || is_fail)
                method_list(k).fails(i) = method_list(k).fails(i)+1;
                continue;
            end            

            %choose the solution with smallest error 
            [index_best, y] = model.choose_best_solution(R1, t1, R, t);
            
            [err errl errt] = model.compute_reprojection(R1, t1, index_best, npt, nlines, XXw, xxn, Xs, Xe, xs, xe);                                    
            
            method_list(k).r(j,i)= y(1);
            method_list(k).t(j,i)= y(2);            
            method_list(k).reproj(j,i) = err;
            method_list(k).reproj_lines(j,i) = errl;
            method_list(k).reproj_pts_lines(j,i) = errt;
            method_list(k).tm1(j, i) = s(1);
%             method_list(k).tm2(j, i) = s(2);
%             method_list(k).tm3(j, i) = s(3);
%             method_list(k).tm4(j, i) = s(4);
        end                    
        showpercent(j,num);    
    end        
    fprintf('\n');
    

    % save result
    for k= 1:length(method_list)
        if (method_list(k).met_type == 1 && nlines < 3)
            continue;
        end
        method_list(k).mean_r(i)= (mean(method_list(k).r(:,i)));
        method_list(k).mean_t(i)= (mean(method_list(k).t(:,i)));
        method_list(k).mean_reproj(i)= (mean(fpix*method_list(k).reproj(:,i)));
        method_list(k).mean_reproj_lines(i)= (mean(fpix*method_list(k).reproj_lines(:,i)));
        method_list(k).mean_reproj_pts_lines(i)= (mean(fpix*method_list(k).reproj_pts_lines(:,i)));
        
        method_list(k).med_r(i)= (median(method_list(k).r(:,i)));
        method_list(k).med_t(i)= (median(method_list(k).t(:,i)));
        method_list(k).med_reproj(i)= median(fpix*method_list(k).reproj(:,i));
        method_list(k).med_reproj_lines(i)= median(fpix*method_list(k).reproj_lines(:,i));
        method_list(k).med_reproj_pts_lines(i)= median(fpix*method_list(k).reproj_pts_lines(:,i));
        method_list(k).avg_t_1(i)= sum(method_list(k).tm1(:,i)) / size(method_list(k).tm1, 1);
%         method_list(k).avg_t_2(i)= sum(method_list(k).tm2(:,i)) / size(method_list(k).tm2, 1);
%         method_list(k).avg_t_3(i)= sum(method_list(k).tm3(:,i)) / size(method_list(k).tm3, 1);
%         method_list(k).avg_t_4(i)= sum(method_list(k).tm4(:,i)) / size(method_list(k).tm4, 1);
    end
end

for k = 1:length(method_list)
    fprintf('%s\n', method_list(k).name);
    method_list(k).fails
end

save('lines_pts_number_timing.mat');

plot_timings;