function [coeffs] = tenscheb_norm2sq(chebyshev_indices)
%**************************************************************************
%
% function [coeffs, chebyshev_indices] = tenscheb_norm2sq(ade, chebyshev_indices)
%
% Calcola i coefficienti di normalizzazione (phi_k, phi_k) per la base
% tensoriale di Chebyshev di grado totale "ade".
%
%**************************************************************************
%
% INPUT:
%
%   ade: 
%           grado totale della base polinomiale
%
% chebyshev_indices: 
%           indici della base tensoriale di Chebyshev. Se chebyshev_indices(s,:) 
%           = [i j k], allora il s-esimo polinomio è: phi_s(x,y,z) 
%           = T_i(x) * T_j(y) * T_k(z) dove T_m(u) è il polinomio di 
%           Chebyshev di grado m
%
%**************************************************************************
%
% OUTPUT:
%
% coeffs: il valore dell'integrale del quadrato di ogni polinomio della base
%         (prodotto scalare del polinomio con se stesso)
%
%**************************************************************************

[~, d] = size(chebyshev_indices);

% Conta in modo vettoriale quanti indici per riga sono diversi da zero
nonzero_counts = sum(chebyshev_indices ~= 0, 2); % N x 1

% Calcola coefficiente pi^d / 2^nonzero_counts
coeffs = (pi^d) ./ (2 .^ nonzero_counts);