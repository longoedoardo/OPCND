function [XI, ETA, T, W_ref, scale_x, num_pts] = reference_prism_quadrature(ade, dbox)
%**************************************************************************
% Genera la griglia di quadratura sul prisma di riferimento (Triangolo x [0,1])
% accoppiando una regola per triangoli basata su Gauss-Jacobi (invece di Dunavant)
% con Gauss-Legendre per l'asse temporale.
%**************************************************************************

    % Numero di punti di Gauss per ogni dimensione lineare
    nGP = ade + 1; 

    % 1. Quadratura sul triangolo di riferimento via Gauss-Jacobi
    % Per assorbire lo Jacobiano (1 - mu1) della mappatura nel triangolo:
    % - mu1 usa Gauss-Jacobi P^(1,0) su [-1, 1]
    % - mu2 usa Gauss-Legendre P^(0,0) su [-1, 1]
    [r1, A1] = gaujac(nGP, 1.0, 0.0);
    [r2, A2] = gaujac(nGP, 0.0, 0.0);

    % Mappatura da [-1, 1] a [0, 1]
    mu1 = 0.5 * r1 + 0.5;   w_mu1 = 0.25 * A1; % 0.5 dal cambio var * 0.5 dallo Jacobiano
    mu2 = 0.5 * r2 + 0.5;   w_mu2 = 0.5  * A2;

    n_mu1 = length(mu1);
    n_mu2 = length(mu2);
    N_spazio = n_mu1 * n_mu2;

    x = zeros(N_spazio, 1);
    y = zeros(N_spazio, 1);
    w = zeros(N_spazio, 1);

    idx = 1;
    for i = 1:n_mu1
        for j = 1:n_mu2
            m1 = mu1(i);
            m2 = mu2(j);

            % Mappatura delle coordinate sul triangolo standard (0,0)-(1,0)-(0,1)
            x(idx) = m1;
            y(idx) = m2 * (1.0 - m1);
            
            % Peso della quadratura spaziale
            w(idx) = w_mu1(i) * w_mu2(j);
            
            idx = idx + 1;
        end
    end

    % 2. Generazione nodi e pesi di Gauss-Legendre per l'asse temporale [0, 1]
    n_p_time = ade + 1; 
    [t, wt] = lgwt_local(n_p_time, 0, 1);

    % 3. Prodotto tensoriale Spazio (Triangolo) x Tempo
    N_tempo = length(wt);
    num_pts = N_spazio * N_tempo;

    XI = zeros(num_pts, 1);
    ETA = zeros(num_pts, 1);
    T = zeros(num_pts, 1);
    W_ref = zeros(num_pts, 1);

    idx = 1;
    for i = 1:N_spazio
        for j = 1:N_tempo
            XI(idx)    = x(i);
            ETA(idx)   = y(i);
            T(idx)     = t(j);
            W_ref(idx) = w(i) * wt(j);
            idx = idx + 1;
        end
    end

    % Fattore di scala fisico per la primitiva lungo X (derivante dalla dbox)
    scale_x = (dbox(2,1) - dbox(1,1)) / 2;
end

%**************************************************************************
%
% GAUSS-LEGENDRE UNIDIMENSIONALE (Golub-Welsch)
%
%**************************************************************************
function [x, w] = lgwt_local(N, a, b)
    N = double(N);
    if N == 1
        x = (a+b)/2; w = b-a;
    else
        i = 1:N-1;
        beta = 0.5 ./ sqrt(1 - (2*i).^(-2));
        Tmat = diag(beta,1) + diag(beta,-1);
        [V, D] = eig(Tmat);
        [x, idx] = sort(diag(D));
        w = 2 * V(1, idx).^2;
    end
    x = x(:); w = w(:);
    x = (a+b)/2 + (b-a)/2 * x;
    w = (b-a)/2 * w;
end

%**************************************************************************
%
% GAUSS-JACOBI UNIDIMENSIONALE (Golub-Welsch)
%
%**************************************************************************
function [x, w] = gaujac(n, alpha, beta)
    if n == 0, x = []; w = []; return; end
    if n == 1
        x = (beta - alpha) / (alpha + beta + 2);
        w = 2^(alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 2);
        return;
    end
    i = (1:n-1)';
    abi = alpha + beta + 2*i;
    aa = zeros(n, 1);
    aa(1) = (beta - alpha) / (alpha + beta + 2);
    aa(2:n) = (beta^2 - alpha^2) ./ (abi .* (abi + 2));
    bb = zeros(n-1, 1);
    bb(1) = 2 * sqrt(1 * (1 + alpha) * (1 + beta) / ((alpha + beta + 2)^2 * (alpha + beta + 3)));
    bb(2:n-1) = 2 ./ (abi(1:end-1) + 2) .* sqrt(i(2:end) .* (i(2:end) + alpha) .* ...
                (i(2:end) + beta) .* (i(2:end) + alpha + beta) ./ ...
                ((abi(1:end-1) + 1) .* (abi(1:end-1) + 3)));
    J = diag(aa) + diag(bb, 1) + diag(bb, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    V = V(:, idx);
    factor = 2^(alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 2);
    w = factor * (V(1, :)').^2;
end