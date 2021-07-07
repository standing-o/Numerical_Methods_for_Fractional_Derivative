function [asn] = numFDeD2(fd2, alpha, a, x, abse, rele)

% fd2 : symbolic second order derivative function | alpha : order | a : start point | 
% x : evaluation point  | abse : absolute error | rele : relative error

g = @(xn) fd2(xn)./((x-xn).^(alpha-1));
asn = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)./gamma(2-alpha);

end