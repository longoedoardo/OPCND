function [vertices,facets, SHP, tri, bbox] = geometryLALOcheap(example)
%**************************************************************************
%
% [vertices,facets, bbox, fv, SHP] = geometryLALO(example)
%
% Genera la geometria di un poliedro 3D di esempio
%
% INPUT:
%   example : numero identificativo del poliedro desiderato
%
% OUTPUT:
%   vertices : Nx3 matrice delle coordinate dei vertici del poliedro
%   facets : Mx3 matrice dei triangoli di superficie (facce)
%   bbox : 1x6 vettore del bounding box [xmin xmax ymin ymax zmin zmax]
%   fv : struct con campi 'vertices' e 'faces', utile per plot
%   SHP : oggetto alphaShape che descrive la superficie del poliedro
%
%**************************************************************************
%
% Viene richiamata la funzione 'provide_vertices' per ottenere i vertici e
% il raggio dell'alphaShape. Viene creato un oggetto alphaShape dai
% vertici, utile per generare le facce del poliedro. In seguito si
% estraggono le facce del poliedro tramite boundaryFacets. Si prepara la
% struttura 'fv' con vertici e facce per plo e funzioni successive. Viene
% infine calcolata la bounding box del poliedro, utile per generare punti
% iniziali per cubatura e poi visualizzazione
%
%**************************************************************************

[vertices, ~, alphashape_radius] = provide_vertices(example);
SHP    = alphaShape(vertices(:,1),vertices(:,2),vertices(:,3),alphashape_radius);
facets = boundaryFacets(SHP);
tri = alphaTriangulation(SHP);

fv.vertices = vertices;
fv.faces = facets;

mins = min(vertices); % Restituisce [minX, minY, minZ]
maxs = max(vertices); % Restituisce [maxX, maxY, maxZ]

% Costruisci il bbox nel formato richiesto [xmin, xmax, ymin, ymax, zmin, zmax]
bbox = [mins(1), maxs(1), mins(2), maxs(2), mins(3), maxs(3)]';
end