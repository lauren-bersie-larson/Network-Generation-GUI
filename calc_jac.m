function [J] = calc_jac(func, x)


% [J] = calc_jac(func, x)
%
% calcs the numerical jacobian given a fnxn by handle -- in netmat
%
% func -- fnxn handle
% x -- input vector
%
% internal parameter for dx (epsilon) affects the calculation significantly
%
% last update -- mon july 23 2012 -- mfh


n = length(x);

fx = feval(func, x);

dx = 1e-8; % larger value than machine eps

xperturb = x;

for i = 1:n

    xperturb(i) = xperturb(i) + dx;   % perturb ith x value

    J(:,i) = (feval(func,xperturb)-fx) / dx; % dF1...n / dxi

    xperturb(i) = x(i); % restore ith xperturb value to original x value

end


end
