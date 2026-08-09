function V_tetr = TensMeshTetrahedron(V, F)

%**************************************************************************
%
% function V_tetr = TensMeshTetrahedron(V, F)
%
% La funzione prende una mesh superficiale definita dai vertici V e dalle 
% facce triangolari F e genera una serie di tetraedri connettendo ogni 
% faccia al baricentro della mesh.
%
% INPUT:
% - V : Matrice (N x 3) dei vertici della superficie.
% - F : Matrice (M x 3) degli indici delle facce (triangoli).
%
% OUTPUT:
% - V_tetr : Tensore (4 x 3 x M). Ogni pagina (:,:,k) rappresenta i 
% 4 vertici del k-esimo tetraedro.
%
% La funzione calcola il baricentro geometrico della nuvola di punti V e lo
% utilizza come vertice comune per tutti i tetraedri. Ogni faccia della 
% superficie diventa la base di un tetraedro il cui apice è il baricentro.
%
%**************************************************************************

% Calcolo del baricentro geometrico della mesh (vettore riga 1x3)
C = mean(V, 1);
numFacce = size(F, 1); % Numero di facce

% Tensore dato da 4 vertici, 3 coordinate (x,y,z), per M tetraedri
V_tetr = zeros(4, 3, numFacce);

% Costruzione dei tetraedri. Per ogni faccia, i primi 3 vertici 
% sono presi da F, il 4o punto è il baricentro
for k = 1:numFacce
    % Estrazione dei vertici della faccia corrente
    V_faccia = V(F(k, :), :);
        
    % Assemblaggio del tetraedro: [Vertice1; Vertice2; Vertice3; Baricentro]
    V_tetr(:, :, k) = [V_faccia; C];
end
end