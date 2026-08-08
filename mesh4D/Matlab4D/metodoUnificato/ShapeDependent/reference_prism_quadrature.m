function [XI_ref, ETA_ref, T_ref, W_ref] = reference_prism_quadrature(ade)
%**************************************************************************
%
% function [XI_ref, ETA_ref, T_ref, W_ref] = reference_prism_quadrature(ade)
%
% Calcola i nodi e i pesi di quadratura di riferimento per il prisma 4D 
% (prodotto cartesiano tra il triangolo 2D della faccia e l'intervallo temporale [0, 1]).
%
%**************************************************************************

% Grado di precisione per la base triangolare (coordinate xi, eta)
nGP_tri = ceil((ade + 2) / 2) + 1;
[tri_points, tri_weights] = TriangleQuadraturePoints(nGP_tri);
n_tri = size(tri_points, 1);

r = tri_points(:, 1);
s = tri_points(:, 2);

% Grado di precisione per la dimensione temporale tau (intervallo [0, 1])
n1D = ceil((ade + 4) / 2);
[tau_nodes, tau_weights] = lgwt(n1D, 0, 1); % Nodi e pesi di Gauss-Legendre su [0, 1]
n_tau = length(tau_nodes);

% Costruzione della griglia come prodotto tensoriale
n_tot = n_tri * n_tau;

XI_ref  = zeros(n_tot, 1);
ETA_ref = zeros(n_tot, 1);
T_ref   = zeros(n_tot, 1);
W_ref   = zeros(n_tot, 1);

idx = 1;
for i = 1:n_tri
    for k = 1:n_tau
        XI_ref(idx)  = r(i);
        ETA_ref(idx) = s(i);
        T_ref(idx)   = tau_nodes(k);   % Vettore colonna singolo per il tempo

        % Il peso totale è il prodotto dei pesi del triangolo e del tempo
        W_ref(idx)   = tri_weights(i) * tau_weights(k);

        idx = idx + 1;
    end
end

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
