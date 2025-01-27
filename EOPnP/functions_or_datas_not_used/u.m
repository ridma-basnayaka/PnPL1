for i = 1 : 50
    test_ordinary_Ximgn;
%     test_planar_Ximgn
    if y(1) > 10
        pause;
    end
    if y(1) > 50
        save error_large;
        pause;
    end
end
