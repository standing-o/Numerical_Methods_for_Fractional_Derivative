function [asn] = RLIeD(f, fd, alpha, a, x, abse, rele)

% f : symbolic function | fd : second order derivative function | alpha : order | a : start point |
% x : evaluation point | abse : absolute error | rele : relative error

temp1 = f(a).*(x-a).^(alpha);
g = @(xn) fd(xn).*(x-xn).^(alpha);
temp2 = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele);
asn = (temp1 + temp2)./gamma(alpha+1);

end