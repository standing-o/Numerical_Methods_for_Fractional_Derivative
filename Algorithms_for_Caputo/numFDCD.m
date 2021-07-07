function numFDCD(f, alpha, a, x, h, abse, rele)

% f : symbolic function | alpha : order | a : start point
% x : evaluation point | h : step size | 
% abse : absolute error | rele : relative error

g = @(xn) imag(f(xn+1*h)./h)./((x-xn).^(alpha));
asn = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)/gamma(1-alpha);

end