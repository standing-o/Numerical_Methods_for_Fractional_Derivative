function [asn] = numFDS(f, alpha, a, x, h, abse, rele)

% f : symbolic function | alpha : order | a : start point
% x : evaluation point | h : step size | 
% abse : absolute error | rele : relative error

g = @(xn)(x-xn).^(-alpha).*((f(xn+h) - f(xn))./h);
asn = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)/gamma(1-alpha);

end
