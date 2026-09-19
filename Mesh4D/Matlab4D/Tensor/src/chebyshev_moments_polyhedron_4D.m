function moments = chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, chebyshev_indices, bbox)

%**************************************************************************
%
% moments = chebyshev_moments_polyhedron_4D(vertici_4D, Hyperfacets, ade, chebyshev_indices, bbox)
%
% Calcola i momenti di Chebyshev su un dominio spazio-temporale 4D usando
% il teorema della divergenza.
%
% Il dominio 4D e' costruito dall'evoluzione lineare di un poliedro 3D tra
% tau = 0 e tau = 1. Le iperfacce laterali sono prismi spazio-temporali
% parametrizzati da due coordinate sulla faccia triangolare e da tau.
%
% Per ogni base di Chebyshev:
%
%   phi(x,y,z,tau) = T_i(x) T_j(y) T_k(z) T_l(tau)
%
% si usa il campo vettoriale:
%
%   F = (Phi_x, 0, 0, 0)
%
% dove Phi_x e' la primitiva di phi rispetto a x. Quindi:
%
%   int_Omega phi dV_4 = int_boundary Phi_x * n_x dS
%
% Poiche' base e tetto hanno normale solo nella direzione tau, non danno
% contributo. Si integrano solo le iperfacce laterali, usando la componente
% locale orientata n_x dS calcolata da PrismQuad4D.
%
% INPUT:
%   vertici_4D        matrice Nv x 4 dei vertici [x, y, z, tau]
%   Hyperfacets       cell array delle iperfacce del dominio 4D
%   ade               grado algebrico di esattezza
%   chebyshev_indices matrice Nm x 4 degli indici della base di Chebyshev
%   bbox              bounding box 4D, formato 2 x 4:
%                     bbox(1,:) = [xmin ymin zmin taumin]
%                     bbox(2,:) = [xmax ymax zmax taumax]
%
% OUTPUT:
%   moments           vettore Nm x 1 dei momenti di Chebyshev sul dominio 4D
%
%**************************************************************************

num_moments = size(chebyshev_indices, 1);
num_facets  = numel(Hyperfacets);

moments = zeros(num_moments, 1);

% Baricentro globale usato da PrismQuad4D per orientare verso l'esterno
baricentro_4D = mean(vertici_4D, 1);

% Regola di quadratura sul prisma di riferimento:
% triangolo 2D x intervallo temporale [0,1].
[XI_ref, ETA_ref, T_ref, W_ref] = reference_prism_quadrature(ade);

% Le prime due iperfacce sono base a tau = 0 e tetto a tau = 1
% Con il campo F = (Phi_x,0,0,0), entrambe hanno n_x = 0

for k = 3:num_facets

    nodes_id = Hyperfacets{k}.Vertices_ID;

    V_prism = vertici_4D(nodes_id, :);

    [XYZTW, WV_X] = PrismQuad4D(V_prism, XI_ref, ETA_ref, T_ref, W_ref, baricentro_4D);

    facet_moments = cubature_tens_chebyshev_facet_4D(XYZTW, WV_X, chebyshev_indices, bbox);

    moments = moments + facet_moments;
end

end