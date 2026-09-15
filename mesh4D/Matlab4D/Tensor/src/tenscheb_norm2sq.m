function coeffs = tenscheb_norm2sq(chebyshev_indices)
%**************************************************************************
%
% function coeffs = tenscheb_norm2sq(chebyshev_indices)
%
% Calcola i coefficienti di normalizzazione (phi_k, phi_k) per la base 
% tensoriale di Chebyshev di grado totale "ade" in 4D.
% Versione interamente vettorializzata (senza cicli for).
%
%**************************************************************************
% Input:
%   chebyshev_indices: Matrice (N x 4) degli esponenti [i, j, k, l]
%
% Output:
%   coeffs: Vettore colonna (N x 1) con il valore dell'integrale T^2
%**************************************************************************

num_indici_non_zero = sum(chebyshev_indices ~= 0, 2);
coeffs = (pi^4) ./ (2 .^ num_indici_non_zero);

end