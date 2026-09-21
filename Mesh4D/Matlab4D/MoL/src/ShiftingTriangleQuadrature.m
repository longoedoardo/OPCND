function [nodes, weights] = ShiftingTriangleQuadrature(V, nodes_ref, weights_ref)

%**************************************************************************
%
% function [nodes, weights] = shiftingTriangleQuadrature( ...
%     V, nodes_ref, weights_ref)
%
% Trasformazione affine dei nodi e dei pesi dal triangolo di riferimento
% al triangolo fisico.
%
% Il triangolo fisico e' definito dai vertici contenuti in V, come 
% V(1,:), V(2,:), V(3,:) La trasformazione utilizzata e':
%
%           x = V1 + xi * (V2-V1) + eta * (V3-V1)
%
%**************************************************************************
%
% INPUT:
%
%   V:
%               Matrice 3 x 3 contenente le coordinate dei vertici del
%               triangolo fisico
%
%   nodes_ref:
%               Matrice N x 2 contenente i nodi di quadratura nel
%               triangolo di riferimento
%
%   weights_ref:
%               Vettore N x 1 contenente i pesi della regola di quadratura
%               sul triangolo di riferimento.
%
%**************************************************************************
%
% OUTPUT:
%
%   nodes:
%               Matrice N x 3 contenente i nodi di quadratura trasformati
%               nel triangolo fisico.
%
%   weights:
%               Vettore N x 1 contenente i pesi della regola di quadratura
%               trasformati nel triangolo fisico.
%
%**************************************************************************

n_points = size(weights_ref, 1);

nodes = zeros(n_points, 3);
weights = zeros(n_points, 1);

% Lati del triangolo reale a partire dal primo vertice.
A = V(2,:) - V(1,:);
B = V(3,:) - V(1,:);

% Trasformazione affine dei punti dal triangolo di riferimento
% al triangolo fisico.
nodes = V(1,:) + nodes_ref(:,1) * A + nodes_ref(:,2) * B;

% Prodotto vettoriale.
cross_product = cross(A, B);

% Area del triangolo fisico.
area = 0.5 * norm(cross_product);

% Trasformazione dei pesi.
weights = area * weights_ref;

end