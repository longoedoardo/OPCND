function moms_facet_raw = cubature_tens_chebyshev_facet_4D(XYZTW, WV_CUB, ade, chebyshev_indices, dbox)
% Calcola l'integrazione polinomiale pura sulla iperfaccia corrente
num_pts = size(XYZTW, 1);
num_indici = size(chebyshev_indices, 1);

box_sum  = dbox(1,:) + dbox(2,:);
box_diff = dbox(2,:) - dbox(1,:);

% Scalatura dei punti sul dominio polinomiale di riferimento [-1, 1]^4
X_ref = (2 * XYZTW - repmat(box_sum, num_pts, 1)) ./ repmat(box_diff, num_pts, 1);

% Valutazione delle basi di Chebyshev 1D
Tx = chebpolys(ade + 1, X_ref(:,1)); 
Ty = chebpolys(ade,     X_ref(:,2));
Tz = chebpolys(ade,     X_ref(:,3));
Tt = chebpolys(ade,     X_ref(:,4));

% Integrazione analitica (Primitiva) lungo la direzione X
scale_x = (dbox(2,1) - dbox(1,1)) / 2;
IntTx = zeros(size(Tx));
IntTx(:, 1) = Tx(:, 2); 
if ade >= 1
    IntTx(:, 2) = 0.25 * Tx(:, 3); 
end
for n = 2:ade
    IntTx(:, n+1) = 0.5 * (Tx(:, n+2)/(n+1) - Tx(:, n)/(n-1));
end
IntTx = IntTx * scale_x; 

% Combinazione tensoriale guidata dagli indici GRLEX
V_poly = zeros(num_pts, num_indici);
for ii = 1:num_indici
    V_poly(:, ii) = IntTx(:, chebyshev_indices(ii, 1)+1) .* ...
                    Ty(:, chebyshev_indices(ii, 2)+1) .* ...
                    Tz(:, chebyshev_indices(ii, 3)+1) .* ...
                    Tt(:, chebyshev_indices(ii, 4)+1);
end

moms_facet_raw = V_poly' * WV_CUB;
end