function [asn] = RLI(f, alpha, a, x, abse, rele)

% f : symbolic function | alpha : order | a : start point |
% x : evaluation point | abse : absolute error | rele : relative error

g = @(xn) f(xn).*((x-xn).^(alpha-1));
asn = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)./gamma(alpha);

end