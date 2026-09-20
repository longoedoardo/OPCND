function [x, w] = lgwt(n, a, b)
%**************************************************************************
%
% function [x, w] = lgwt(n, a, b)
% 
% La funzione calcola i nodi (x) e i pesi (w) della quadratura di 
% Gauss-Legendre nell'intervallo generico [a, b]. Corrisponde al calcolo di
% nodi e pesi della quadratura di Gauss-Jacobi con alpha = beta = 0.
%
%**************************************************************************

i = 1:n-1;
a_coef = ones(n,1);
b_coef = i ./ sqrt(4*i.^2 - 1);

CM = diag(b_coef, 1) + diag(b_coef, -1);
[V, L] = eig(CM);

[x, ind] = sort(diag(L));
V = V(:, ind);

w = 2 * (V(1, :)').^2;

x = (a*(1 - x) + b*(1 + x)) / 2;
w = w * (b - a) / 2;
end