function moments = chebyshev_moments_polyhedron_4D(vertici_4D, facets, ...
                                      ade, chebyshev_indices, bbox, method)

%**************************************************************************
%
% moments = chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, 
%                                           chebyshev_indices, bbox, method)
%
% Calcola i momenti di Chebyshev su un dominio spazio-temporale 4D usando
% il teorema della divergenza.
%
% Per ogni phi(x,y,z,tau) = T_i(x) T_j(y) T_k(z) T_l(tau) si usa il campo 
% vettoriale F = (Phi_x, 0, 0, 0) dove Phi_x e' la primitiva di phi 
% rispetto a x. Poiche' base e tetto hanno normale solo nella direzione 
% tau, non danno contributo. Si integrano solo le iperfacce laterali, 
% usando la componente locale orientata n_x dS.
%
%**************************************************************************
%
% INPUT:
%
%   vertici_4D:        
%               Matrice Nv x 4 dei vertici [x, y, z, tau]
%
%   facets:
%               Matrice Nf x 3 contenente gli indici delle facce triangolari
%               della mesh tridimensionale. Ogni faccia individua
%               un'iperfaccia laterale del dominio spazio-temporale
%
%   ade:
%               Grado massimo della base polinomiale di Chebyshev per la
%               quale vengono calcolati i momenti
%
%   chebyshev_indices: 
%               Matrice Nm x 4 degli indici della base di Chebyshev
%
%   bbox:              
%               bounding box 4D, formato 2 x 4:
%
%                     bbox(1,:) = [xmin ymin zmin taumin]
%                     bbox(2,:) = [xmax ymax zmax taumax]
%
%   method:
%               Metodo di quadratura sul prisma di riferimento:
%               - 'DCC'  = Dunavant--Clenshaw--Curtis
%               - 'DGL'  = Dunavant--Gauss--Legendre
%               - 'GJCC' = Gauss--Jacobi--Clenshaw--Curtis
%               - 'GJL'  = Gauss--Jacobi--Gauss--Legendre
%
%**************************************************************************
%
% OUTPUT:
%   moments:           
%           Vettore dei momenti di Chebyshev sul dominio 4D
%
%**************************************************************************

num_moments = size(chebyshev_indices, 1);
num_facets  = size(facets, 1);
num_vertici = size(vertici_4D, 1) / 2;

moments = zeros(num_moments, 1);

% Costruzione della regola di quadratura sul prisma di riferimento
[XI_ref, ETA_ref, T_ref, W_ref] = PrismQuadrature(ade, method);

for k = 1:num_facets

    n1 = facets(k,1);
    n2 = facets(k,2);
    n3 = facets(k,3);

    n1_new = n1 + num_vertici;
    n2_new = n2 + num_vertici;
    n3_new = n3 + num_vertici;

    V_prism = vertici_4D([n1 n2 n3 n1_new n2_new n3_new], :);

    [XYZTW, WV_X] = ShiftingPrismQuadrature(V_prism, XI_ref, ETA_ref, T_ref, W_ref);
        
    % Se la componente x della normale e' nulla,
    % l'iperfaccia non da' alcun contributo.
    if max(abs(WV_X)) < 1e-14
        continue
    end

    facet_moments = cubature_tens_chebyshev_facet_4D(XYZTW, WV_X, chebyshev_indices, bbox);

    moments = moments + facet_moments;
end

end