function [x_opt,fval] = minimize_kernel_diff_L2(psfr,kernel_fun,x0)

    % Minimize the difference between a given Point-spread-function and a functional guess of the kernel.
    % Return the optimal paramters and the L2 norm at the found optimum.
    
    diff = @(x) sum((psfr - kernel_fun(x)).^2,'all');
    [x_opt,fval] = fminsearch(diff,x0);
end