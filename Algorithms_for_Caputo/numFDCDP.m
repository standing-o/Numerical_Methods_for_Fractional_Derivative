function [asn] = numFDCDP(f, alpha, a, x, h1, h2, abse, rele)

% f : symbolic function | alpha : order | a : start point | x : evaluation point 
% h1 : stepsize for first order | h2 : stepsize for second order
% abse : absolute error | rele : relative error

temp1 = imag(f(a+1j.*h1)./h1).*(x-a).^(1-alpha)./gamma(2-alpha);
I = h2./sqrt(2).*(1+1j);
g = @(xn) imag((f(xn+I) + f(xn-I))./((h2)^2)).*((x-xn).^(1-alpha));
temp2 = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)/gamma(2-alpha);

asn = temp1 + temp2;

end
 