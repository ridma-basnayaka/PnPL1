clear; 
clc;
close all;
warning off;
rng('default');
prepare_paths();

% experimental parameters

%focal length in pixels
fpix = 800;
%varying noise level
nls = [0.5:0.5:2];
%length of a line segment in 3D 
segLenBig = 3;
%scale of shift between line segments projected and present in the model
segLenShiftScale = 0.5*segLenBig;

%number of random trials per iteration
num = 100;
%number of lines
nlines = 10;
%number of points
npt = 0;

%methods chosen for comparison   
% method_names = { 'RPnL', 'DLT', 'EPnPL', 'OPnPL', 'EOPnPL','EOPnP*'};
% method_names = { 'Mirzaei','ASPnL','LPnL\_Bar\_ENull','Ansar','RPnL', 'SRPnL','DLT', 'EPnPL', 'OPnPL', 'EOPnPL'};
% method_names = {'EPnPL','OPnPL','EOPnPL','EOPnPL-BJ-0','EOPnPL-BJ-1','EOPnPL-BJ-2','EOPnPL-BJ-3','EOPnPL-BJ-4'};
method_names = {'EPnPL','OPnPL','EOPnPL','EOPnPL-BJ-0'};


is_planar = 0;
method_list = model.setup_methods(method_names, num, length(nls), is_planar);


%experiments
for i= 1:length(nls)
            
    fprintf('nls = %d: ',nls(i));
%noise in pix, added to point projections, std. dev.
    nll = nls(i);
%noise in pix, added to line segment endpoints' projections, std. dev.
    nlp = nls(i);
    
    index_fail = [];
    
    for k = 1:length(method_list)
        method_list(k).fails(i) = 0;
    end   

    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;       
        
        [R, t, XLTw, XXw, xxn, xs,xe,Xs,Xe] = model.setup_3d_scene(npt, nlines, segLenBig, segLenShiftScale, nll, nlp, fpix);
                             
            

        for k = 1:length(method_list)
            
            % pose estimation
            [R1, t1, is_fail, s] = model.evaluate_method(method_list(k), Xs, Xe, xs, xe, XXw, xxn, XLTw); 
            if (size(R1,1) < 3 || is_fail)
                method_list(k).fails(i) = method_list(k).fails(i)+1;
                continue;
            end            

            %choose the solution with smallest error 
            [index_best, y] = model.choose_best_solution(R1, t1, R, t);
            
            if (y(2) > 100)
                y;
            end
            [err errl errt] = model.compute_reprojection(R1, t1, index_best, npt, nlines, XXw, xxn, Xs, Xe, xs, xe);                                    
            
            method_list(k).r(j,i)= y(1);
            method_list(k).t(j,i)= y(2);            
            method_list(k).reproj(j,i) = err;
            method_list(k).reproj_lines(j,i) = errl;
            method_list(k).reproj_pts_lines(j,i) = errt;
            method_list(k).tm(j, i) = s;
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
        method_list(k).avg_t(i)= sum(method_list(k).tm(:,i)) / size(method_list(k).tm, 1);
    end
end

for k = 1:length(method_list)
    fprintf('%s\n', method_list(k).name);
    method_list(k).fails
end

save('e_lines_noise.mat');

e_plot_lines_noise;