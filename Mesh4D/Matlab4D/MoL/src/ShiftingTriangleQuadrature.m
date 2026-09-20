function [nodes, weights] = ShiftingTriangleQuadrature(V, nodes_ref, weights_ref)

%******************************************************************************
%
% function [nodes, weights] = shiftingTriangleQuadrature( ...
%     V, nodes_ref, weights_ref)
%
% Trasformazione affine dei nodi e dei pesi dal triangolo di riferimento
% al triangolo fisico.
%
% Il triangolo fisico e' definito dai vertici contenuti in V:
%
%       V(1,:), V(2,:), V(3,:)
%
% La trasformazione utilizzata e':
%
% x = V1 + xi * (V2-V1) + eta * (V3-V1)
%
%******************************************************************************

n_points = size(weights_ref, 1);

nodes = zeros(n_points, 3);
weights = zeros(n_points, 1);

% Lati del triangolo reale a partire dal primo vertice.
A = V(2,:) - V(1,:);
B = V(3,:) - V(1,:);

% Trasformazione affine dei punti dal triangolo di riferimento
% al triangolo fisico.
for i = 1:n_points

    nodes(i,:) = V(1,:) ...
        + nodes_ref(i,1) * A ...
        + nodes_ref(i,2) * B;

end

% Prodotto vettoriale.
cross_product = cross(A, B);

% Area del triangolo fisico.
area = 0.5 * norm(cross_product);

% Trasformazione dei pesi.
weights = area * weights_ref;

end