function [x_opt,fval] = minimize_kernel_diff_L2(psfr,kernel_fun,x0)
diff = @(x) sum((psfr - kernel_fun(x)).^2,'all');
[x_opt,fval] = fminsearch(diff,x0);
end