function [asn] = numFDCD2(f, alpha, a, x, h, abse, rele)

% f : symbolic function | alpha : order | a : start point | x : evaluation point 
% h : stepsize | abse : absolute error | rele : relative error

I = h*1/sqrt(2)*(1+1j);
g = @(xn) (x-xn).^(-(alpha-1)).*imag((f(xn+I)+f(xn-I))/((h)^2));
asn = integral(g, a, x, 'AbsTol', abse, 'RelTol', rele)/gamma(2-alpha);

end