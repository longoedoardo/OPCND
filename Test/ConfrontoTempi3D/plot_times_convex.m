clc
clear
close all

projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'Matlab3D/src'));
addpath(fullfile(projectRoot, 'TetraMethod'));
addpath(fullfile(projectRoot, 'examples'));

% =========================================================================
% Parametri e caricamento mesh
% =========================================================================

vertices = load('convex_vertex.dat');
facets   = load('convex_tri.dat');

ade_values = 1:19;

time_GJ = zeros(size(ade_values));
time_D  = zeros(size(ade_values));

% =========================================================================
% Benchmark costruzione delle regole di quadratura
% =========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('       BENCHMARK COSTRUZIONE REGOLE DI QUADRATURA\n');
fprintf('============================================================\n');

for k = 1:length(ade_values)

    ade = ade_values(k);

    fprintf('\nADE = %d\n', ade);

    % ---------------------------------------------------------------------
    % Gauss-Jacobi
    % ---------------------------------------------------------------------

    time_GJ(k) = timeit(@() OPC3D(ade, vertices, facets, "GJ"));

    fprintf('  Gauss-Jacobi : %.6e s\n', time_GJ(k));

    % ---------------------------------------------------------------------
    % Dunavant
    % ---------------------------------------------------------------------

    time_D(k) = timeit(@() OPC3D(ade, vertices, facets, "D"));

    fprintf('  Dunavant     : %.6e s\n', time_D(k));

end

% =========================================================================
% Plot dei tempi
% =========================================================================

figure('Color', 'w');

plot(ade_values, time_GJ, '-o', ...
     'LineWidth', 1.5, ...
     'MarkerSize', 5);

hold on;

plot(ade_values, time_D, '-s', ...
     'LineWidth', 1.5, ...
     'MarkerSize', 5);

hold off;

grid on;
box on;

xlabel('$ade$', 'Interpreter', 'latex');
ylabel('Tempo [s]', 'Interpreter', 'latex');

legend({'Gauss-Jacobi', 'Dunavant'}, ...
       'Location', 'northwest', ...
       'Interpreter', 'latex');

xticks(ade_values);

set(gca, ...
    'FontSize', 11, ...
    'TickLabelInterpreter', 'latex');

ytickformat('%.1e');

xlim([1 19]);