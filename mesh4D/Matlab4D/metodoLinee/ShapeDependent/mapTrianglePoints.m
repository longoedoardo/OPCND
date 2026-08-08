function [nodifisici, pesifisici, norm_ext] = mapTrianglePoints(nodi_rif, pesi_rif, V)
%**************************************************************************
%
% function [nodi_fisici, pesi_fisici, norm_ext] = mapTrianglePoints(nodi_rif, pesi_rif, V)
%
% Esegue una trasformazione affine per trasportare i punti di integrazione
% definiti sul triangolo di riferimento (vertici (0,0), (1,0), (0,1)) su un
% triangolo generico in 3D definito dai vertici in V, e calcola la normale
% esterna al triangolo tramite prodotto vettoriale dei due lati.
%
% INPUT:
%     nodi_rif  : Matrice (N x 2) dei nodi nel dominio di riferimento (r,s).
%     pesi_rif  : Vettore colonna (N x 1) dei pesi di quadratura associati.
%                   Si assume che sum(pesi_rif) = 1/2 (area del rif.).
%     V         : Matrice (3 x 3) dei vertici del triangolo target.
%                   V(1,:) è l'origine locale (A); V(2,:), V(3,:) sono (B, C).
%
% OUTPUT:
%     nodi_fisici : Matrice (N x 3) dei nodi mappati nello spazio fisico (x, y, z).
%     pesi_fisici : Vettore (N x 1) dei pesi riscalati per l'area del target.
%     norm_ext    : Versore (1 x 3) normale al piano del triangolo, con verso
%                   determinato dalla regola della mano destra applicata
%                   all'ordine dei vertici V(1,:)->V(2,:)->V(3,:). Per un
%                   poliedro con facce orientate in senso antiorario viste
%                   dall'esterno (convenzione richiesta da
%                   chebyshev_moments_polyhedron), questa è la normale
%                   uscente.
%
%**************************************************************************

% Vettori spigolo dal primo vertice
v1 = V(2,:) - V(1,:);
v2 = V(3,:) - V(1,:);

% Prodotto vettoriale: modulo = 2 * Area fisica, direzione = normale
nvec = cross(v1, v2);
area2 = norm(nvec);

if area2 < eps
    warning('mapTrianglePoints: triangolo degenere o quasi-piatto.');
    norm_ext = [0 0 0];
else
    norm_ext = nvec / area2;
end

% Mappa affine: P(r,s) = V1 + r*v1 + s*v2
r = nodi_rif(:,1);
s = nodi_rif(:,2);
nodifisici = V(1,:) + r*v1 + s*v2;

% I pesi scalano con il rapporto (Area_fisica / Area_riferimento).
% Area_fisica = area2/2 ; Area_riferimento = 1/2  =>  fattore = area2
pesifisici = pesi_rif * area2;

end