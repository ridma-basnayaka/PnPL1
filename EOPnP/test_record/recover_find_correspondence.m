syms A B C s1 s2 s3 real
fs1 = (B^2+C^2+A^2) * s1^2 + 4*C*s1 + C^2;
fs2 = (B^2+C^2+A^2) * s2^2 - 4*B*s2 + B^2;
fs3 = (B^2+C^2+A^2) * s3^2 + 4*A*s3 + A^2;
sol_s1 = solve(fs1);
sol_s2 = solve(fs2);
sol_s3 = solve(fs3);

flag = 1;   %标志位，用于跳出多重循环
for i = 1:2 
    s1 = sol_s1(i);
    for j = 1 : 2
        s2 = sol_s2(j);
        for k = 1 : 2
            s3 = sol_s3(k);
            f1 = (4 * s3) + A * (1 + s1 * s1 + s2 * s2 + s3 * s3);
            f2 = (4 * s2) - B * (1 + s1 * s1 + s2 * s2 + s3 * s3);
            f3 = (4 * s1) + C * (1 + s1 * s1 + s2 * s2 + s3 * s3);
            if(simplify(f1 * f1 + f2 * f2 + f3 * f3) == 0)
                disp('符合条件的一组s1,s2,s3: ')
                ccode(s1)
                ccode(s2)
                ccode(s3)
            end
        end
    end
end