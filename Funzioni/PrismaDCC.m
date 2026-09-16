function [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDCC(ade)

%**************************************************************************
%
% [XI_ref, ETA_ref, T_ref, W_ref] = PrismaDCC(ade)
%
% Costruisce la regola di quadratura tensoriale Dunavant--Clenshaw--Curtis
% (DCC) sul prisma triangolare di riferimento
%
%       P_ref = { (xi,eta,tau) : xi>=0, eta>=0, xi+eta<=1, 0<=tau<=1 },
%
% ottenuta come prodotto tensoriale tra:
%
%   - la regola simmetrica di Dunavant sul triangolo di riferimento 
%
%   - la regola di Clenshaw--Curtis sull'intervallo temporale [0,1]
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

rule_tri = ade;
n_tri = dunavant_order_num(rule_tri);
[xyrif, wrif] = dunavant_rule(rule_tri, n_tri);

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

n1D = ade+1;
[tau_nodes, tau_weights] = clenshaw_curtis(n1D, 0, 1);



n_tau = numel(tau_nodes);

XI_ref  = kron(ones(n_tau, 1), r);
ETA_ref = kron(ones(n_tau, 1), s);
T_ref   = kron(tau_nodes, ones(n_tri, 1));
W_ref   = kron(tau_weights, tri_weights);

end