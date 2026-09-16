function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDGL(ade)

%**************************************************************************
%
% function [XI_ref, ETA_ref, T_ref, W_ref] = DGL(ade)
%
% Costruisce la regola di quadratura tensoriale Dunavant--Gauss--Legendre
% (DGL) sul prisma triangolare di riferimento
%
%       P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
%
% ottenuta come prodotto tensoriale tra:
%
%   - una regola simmetrica di Dunavant sul triangolo di riferimento
%     (funzioni dunavant_order_num e dunavant_rule), che sfrutta le
%     simmetrie del gruppo diedrale D3;
%
%   - la stessa regola di Gauss--Legendre sull'intervallo [0,1] impiegata
%     in GJL (funzione lgwt).
%
%**************************************************************************
%
%   INPUT:
%       ade      : (scalare intero) grado di esattezza algebrica target che si
%                   vuole raggiungere sul prisma.
%
%   OUTPUT:
%       XI_ref, ETA_ref  : (Nq x 1) coordinate baricentriche (xi,eta) dei nodi
%                           di quadratura sul triangolo di riferimento,
%                           ripetute per ciascun nodo temporale.
%       T_ref            : (Nq x 1) coordinata temporale tau in [0,1] di
%                           ciascun nodo.
%       W_ref            : (Nq x 1) peso complessivo di ciascun nodo sul
%                           prisma di riferimento, dato dal prodotto tra il
%                           peso spaziale e il peso temporale.
%
%**************************************************************************

%**************************************************************************
%
% Inizio Quadratura Spaziale
%
%**************************************************************************

% dunavant_order_num restituisce il numero totale di punti fisici n_tri
% che la regola di grado rule_tri produrrà una volta espansa (somma delle
% molteplicità 1, 3 o 6 di ciascuna orbita di simmetria).
n_tri = dunavant_order_num(ade);

% dunavant_rule espande i subordini tabulati nell'elenco completo dei
% nodi e dei pesi, applicando le simmetrie del gruppo D3 (permutazioni
% cicliche e riflessioni delle coordinate baricentriche). Restituisce le
% coordinate come matrice 2 x n_tri e i pesi come vettore di lunghezza n_tri.
[xyrif, wrif] = dunavant_rule(ade, n_tri);

% Portiamo i nodi spaziali in formato colonna
r = xyrif(1, :).';
s = xyrif(2, :).';

if isrow(wrif)
    wrif = wrif.';
end

tri_weights = 0.5 * wrif;

%**************************************************************************
%
% Inizio Quadratura Temporale
%
%**************************************************************************

n1D = ceil((ade + 1) / 2);
[tau_nodes, tau_weights] = lgwt(n1D, 0, 1);

n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end