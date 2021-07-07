function [asn] = RLICD(f, alpha, a, x, h, abse, rele)

% f : symbolic function | alpha : order | a : start point | h : stepsize
% x : evaluation point | abse : absolute error | rele : relative error

temp1 = f(a).*(x-a).^(alpha);
g = @(xn) imag(f(xn+1j*h)./h).*(x-xn).^(alpha);
temp2 = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele);
asn = (temp1 + temp2)./gamma(alpha+1);

end