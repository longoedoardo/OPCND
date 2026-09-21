function [tau_nodi, tau_pesi] = ClenshawCurtisTime(n_tau)

%**************************************************************************%
% function [tau_nodi, tau_pesi] = ClenshawCurtisTime(n_tau)
%
% Costruisce i nodi di Chebyshev-Lobatto e i relativi pesi della formula
% di quadratura di Clenshaw-Curtis sull'intervallo temporale [0,1].
% I nodi sono ottenuti mediante una trasformazione affine dei nodi di
% Chebyshev-Lobatto definiti sull'intervallo [-1,1].
%
% INPUT:
%     n_tau     : Numero di nodi di quadratura temporale.
%
% OUTPUT:
%     tau_nodi  : Vettore colonna (n_tau x 1) dei nodi di Chebyshev-Lobatto
%                 nell'intervallo [0,1], ordinati in senso crescente.
%                 I due estremi 0 e 1 sono inclusi.
%     tau_pesi  : Vettore colonna (n_tau x 1) dei pesi della formula di
%                 quadratura di Clenshaw-Curtis associati ai nodi tau_nodi.
%                 I pesi sono normalizzati sull'intervallo [0,1], pertanto
%                 sum(tau_pesi) = 1.
%
%**************************************************************************

% La formula di Clenshaw-Curtis utilizza i nodi di Chebyshev-Lobatto
% x_i = cos(i*pi/N) con i = 0,...,N, sull'intervallo [-1,1]. 
% Se n_tau è il numero totale di nodi, allora il numero di sottointervalli 
% è N = n_tau - 1

num_intervalli = n_tau - 1;
theta = pi * (0:num_intervalli)' / num_intervalli;
nodiChebLob = cos(theta);

% Trasformazione affine da [-1,1] a [0,1] posto x -> (x+1)/2
% I nodi vengono successivamente riordinati in ordine crescente
tau_nodi = flipud((nodiChebLob + 1) / 2);

% Nella formula dei pesi compaiono i termini pari, cos(2*j*theta)
% Sono quindi necessari gli indici j = 1,...,floor(N/2).
indici_serie = (1:floor(num_intervalli/2))';

% Coefficienti che moltiplicano i termini della serie
coefficienti_serie = 2 * ones(size(indici_serie));

% Se N è pari, l'ultimo termine della serie corrisponde a j = N/2
% In questo caso il suo coefficiente deve essere 1 anziché 2
if mod(num_intervalli, 2) == 0
    coefficienti_serie(end) = 1;
end

% Costruiamo una matrice in cui ogni riga corrisponde a un nodo theta_i
% e ogni colonna corrisponde a un indice j. L'elemento (i,j) è
% [cos(2*j*theta_i)] / [1 - 4*j^2] che compare nella formula dei
% pesi di Clenshaw-Curtis

termini_serie = ...
    cos(2 * theta * indici_serie') ...
    ./ (1 - 4 * indici_serie'.^2);

% Sommiamo i termini della serie pesandoli con i coefficienti
% precedentemente definiti
somma_serie = termini_serie * coefficienti_serie;

% Nella formula dei pesi, i nodi agli estremi theta = 0 e theta = pi
% hanno un fattore 1/2. Per tutti gli altri nodi il fattore vale 1
fattore_estremi = ones(n_tau, 1);
fattore_estremi([1, end]) = 0.5;

% Formula dei pesi sull'intervallo [-1,1]:
%
% w_i = 1/N * fattore_i * (1 + somma_serie_i).
%
% Poiché i nodi sono stati trasformati nell'intervallo [0,1],
% la corrispondente normalizzazione viene mantenuta nella formula

tau_pesi = ...
    (1 / num_intervalli) ...
    * fattore_estremi ...
    .* (1 + somma_serie);

% I nodi di Chebyshev erano inizialmente ordinati da 1 a 0
% Dopo flipud, i nodi tau_nodi sono invece ordinati da 0 a 1
% I pesi devono avere lo stesso ordinamento dei nodi
tau_pesi = flipud(tau_pesi);

end