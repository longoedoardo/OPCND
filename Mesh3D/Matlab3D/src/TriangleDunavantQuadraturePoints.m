function [nodi_rif, pesi_rif] = TriangleDunavantQuadraturePoints(rule)
%********************************************************************************
%
% function [nodi_rif, pesi_rif] = TriangleDunavantQuadraturePoints(rule)
%
% Calcolo dei nodi e dei pesi di una regola di quadratura di Dunavant
% sul triangolo di riferimento.
%
% INPUT:
%   rule       - (int) Grado della regola di quadratura Dunavant.
%
% OUTPUT:
%   nodi_rif   - (Np x 2) Coordinate dei nodi di quadratura sul triangolo
%                di riferimento.
%   pesi_rif   - (Np x 1) Pesi di quadratura sul triangolo di riferimento.
%
% Per poter utilizzare la seguente funzione, sono necessari gli script
% contenuti nella cartella "Dunavant".
%
%********************************************************************************

thisFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(thisFolder, 'Dunavant'));

% Calcolo indice della regola
order_num = dunavant_order_num(rule);

% Calcolo coordinate e pesi sul triangolo di riferimento
[xyrif, wrif] = dunavant_rule(rule, order_num);

% Formatto le coordinate come matrice Np x 2
nodi_rif = xyrif.';

% Formatto i pesi come vettore colonna
if isrow(wrif)
    wrif = wrif.';
end

pesi_rif = wrif;

end