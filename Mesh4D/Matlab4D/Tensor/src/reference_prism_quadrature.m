function [XI_ref, ETA_ref, T_ref, W_ref] = reference_prism_quadrature(ade)

%**************************************************************************
%
% [XI_ref, ETA_ref, T_ref, W_ref] = reference_prism_quadrature(ade)
%
% Costruisce la regola di quadratura di riferimento sul prisma
% triangolare spazio-temporale, ottenuta come prodotto tensoriale tra
% una regola sul triangolo di riferimento e una regola di Gauss-Legendre
% sull'intervallo temporale [0,1].
%
%**************************************************************************

%**************************************************************************
%
% Quadratura sul triangolo di riferimento
%
%**************************************************************************

nGP_tri = ceil((ade + 2) / 2);

[tri_points, tri_weights] = TriangleQuadraturePoints(nGP_tri);

r = tri_points(:, 1);
s = tri_points(:, 2);

%**************************************************************************
%
% Quadratura di Gauss-Legendre nella direzione temporale
%
%**************************************************************************

n1D = ceil((ade + 4) / 2);

[tau_nodes, tau_weights] = lgwt(n1D, 0, 1);

%**************************************************************************
%
% Prodotto tensoriale delle due regole
%
%**************************************************************************

n_tri = numel(r);
n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end


function [x, w] = lgwt(n, a, b)
%**************************************************************************
% Calcola i nodi (x) e i pesi (w) di Gauss-Legendre nell'intervallo [a, b].
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
