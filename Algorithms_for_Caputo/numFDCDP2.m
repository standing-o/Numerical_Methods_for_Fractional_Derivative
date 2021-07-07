function [asn] = numFDCDP2(f, alpha, a, x, h1, h2, h3, abse, rele)

% f : symbolic function | alpha : order | a : start point | x : evaluation point 
% h1, h2, h3 : stepsize | abse : absolute error | rele : relative error

I = h2*1/sqrt(2)*(1+1j);
temp1 = (imag((f(a+I)+f(a-I))/((h2)^2)).*(x-a).^(2-alpha))./gamma(3-alpha);
g = @(xn) (3.*(imag(f(xn-1j.*h3) - f(xn+1j.*h3)) + 2.*h3.*imag(f(xn+1j.*h1)./h1)).*((x-xn).^(2-alpha)))/(h3.^3);
temp2 = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)./gamma(3-alpha);
asn = temp1 + temp2;

end