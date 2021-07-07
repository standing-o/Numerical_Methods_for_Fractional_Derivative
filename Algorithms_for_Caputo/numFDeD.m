function [asn] = numFDeD(fd, alpha, a, x, abse, rele)

% df : symbolic derivative function | alpha : order | a : start point | 
% x : evaluation point |  abse : absolute error | rele : relative error

g = @(xn)(x-xn).^(-alpha).*fd(xn);
asn = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)/gamma(1-alpha);

end