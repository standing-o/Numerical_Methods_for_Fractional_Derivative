function [asn] = numFDeDP(fd1, fd2, alpha, a, x, abse, rele)

% fd1 : symbolic first order derivative function | fd2 : symbolic second order derivative function | 
% alpha : order | a : start point | x : evaluation point 
% abse : absolute error | rele : relative error

temp1 = fd1(a).*(x-a).^(1-alpha)./gamma(2-alpha);
g = @(xn) fd2(xn).*((x-xn).^(1-alpha));
temp2 = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)/gamma(2-alpha);
asn = temp1 + temp2;

end