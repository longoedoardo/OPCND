function [x, w] = clenshaw_curtis(n, a, b)
%**************************************************************************
%
% function [x, w] = clenshaw_curtis(n, a, b)
% 
% La funzione calcola i nodi (x) e i pesi (w) della quadratura di 
% Clenshaw-Curtis nell'intervallo generico [a, b].
%
% I nodi sono ottenuti come trasformazione affine dei nodi
%
%       x_j = cos(j*pi/(n-1)),
%
% definiti sull'intervallo [-1,1].
%
%**************************************************************************

if n == 1
    x = (a + b) / 2;
    w = b - a;
    return;
end

% Creazione dei nodi di Chebyshev-Lobatto

j = (0:n-1)';
x = cos(pi * j / (n-1));

% Calcolo dei pesi di Clenshaw--Curtis
w = zeros(n, 1);

N = n - 1;
j_int = (1:N-1)';

if mod(N, 2) == 0

    % Caso N pari
    w(1) = 1 / (N^2 - 1);
    w(end) = w(1);

    w_int = ones(N-1, 1);

    for k = 1:(N/2 - 1)
        w_int = w_int ...
            - 2 * cos(2*pi*k*j_int/N) ...
            / (4*k^2 - 1);
    end

    w_int = w_int ...
        - (-1).^j_int / (N^2 - 1);

    w(2:N) = (2/N) * w_int;

else

    % Caso N dispari
    w(1) = 1 / N^2;
    w(end) = w(1);

    w_int = ones(N-1, 1);

    for k = 1:((N-1)/2)
        w_int = w_int ...
            - 2 * cos(2*pi*k*j_int/N) ...
            / (4*k^2 - 1);
    end

    w(2:N) = (2/N) * w_int;

end

% Mettiamo in ordine crescente i nodi
[x, ind] = sort(x);
w = w(ind);

% Applicata una trasformazione dall'intervallo [-1,1] all'intervallo [a,b]
% per ricondursi al caso generale
x = (a + b)/2 + (b - a)/2 * x;
w = w * (b - a)/2;

end