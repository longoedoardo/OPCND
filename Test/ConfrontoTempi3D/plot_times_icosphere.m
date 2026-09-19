
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'Matlab3D/src'));
addpath(fullfile(projectRoot, 'TetraMethod'));
addpath(fullfile(projectRoot, 'examples'));

%**************************************************************************
% Parametri e caricamento mesh
%**************************************************************************

vertices = load('icosphere_vertex.dat');

facets = load('icosphere_tri.dat');

ade_values = 1:19;

time_GJ = zeros(size(ade_values));
time_D  = zeros(size(ade_values));

%**************************************************************************
% Benchmark costruzione delle regole di quadratura
%**************************************************************************

fprintf('\n');
fprintf('============================================================\n');
fprintf('       BENCHMARK COSTRUZIONE REGOLE DI QUADRATURA\n');
fprintf('============================================================\n');

for k = 1:length(ade_values)

    ade = ade_values(k);

    fprintf('\nADE = %d\n', ade);

    %**********************************************************************
    % Gauss-Jacobi
    %**********************************************************************

    time_GJ(k) = timeit(@() OPC3D(ade, vertices, facets, "GJ"));

    fprintf('  Gauss-Jacobi : %.6e s\n', time_GJ(k));

    %**********************************************************************
    % Dunavant
    %**********************************************************************

    time_D(k) = timeit(@() OPC3D(ade, vertices, facets, "D"));

    fprintf('  Dunavant     : %.6e s\n', time_D(k));

end




%**************************************************************************
% Plot dei tempi
%**************************************************************************
figure('Color', 'w');

%**************************************************************************
% Plot principale
%**************************************************************************
ax_main = axes;

plot(ax_main, ade_values, time_GJ, '-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5);
hold(ax_main, 'on');

plot(ax_main, ade_values, time_D, '-s', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5);

grid(ax_main, 'on');
box(ax_main, 'on');

xlabel(ax_main, '$ade$', 'Interpreter', 'latex');
ylabel(ax_main, 'Tempo [s]', 'Interpreter', 'latex');

legend(ax_main, {'Gauss-Jacobi', 'Dunavant'}, ...
    'Location', 'northwest', ...
    'Interpreter', 'latex');

xticks(ax_main, ade_values);

set(ax_main, ...
    'FontSize', 11, ...
    'TickLabelInterpreter', 'latex');

ytickformat(ax_main, '%.1e');
xlim(ax_main, [1 19]);

%**************************************************************************
% Riquadro della zona ingrandita: ade = 1,...,7
%**************************************************************************
idx = (ade_values >= 1) & (ade_values <= 7);

y_zoom = [time_GJ(idx), time_D(idx)];
ymin = min(y_zoom);
ymax = max(y_zoom);

y_margin = 0.03 * (ymax - ymin);

rectangle(ax_main, ...
    'Position', [1, ymin-y_margin, 6, ...
                 (ymax-ymin)+2*y_margin], ...
    'EdgeColor', 'k', ...
    'LineWidth', 0.5, ...
    'LineStyle', '-');

%**************************************************************************
% Inset
%**************************************************************************
ax_inset = axes('Position', [0.20 0.28 0.30 0.25]);

plot(ax_inset, ade_values(idx), time_GJ(idx), '-o', ...
    'LineWidth', 1.4, ...
    'MarkerSize', 4);

hold(ax_inset, 'on');

plot(ax_inset, ade_values(idx), time_D(idx), '-s', ...
    'LineWidth', 1.4, ...
    'MarkerSize', 4);

grid(ax_inset, 'on');
box(ax_inset, 'on');

xlim(ax_inset, [1 7]);
xticks(ax_inset, 1:7);

ylim(ax_inset, [ymin ymax]);

set(ax_inset, ...
    'FontSize', 8, ...
    'TickLabelInterpreter', 'latex');

ytickformat(ax_inset, '%.2e');


%**************************************************************************
% Collegamento tra il riquadro e l'inset
%**************************************************************************
drawnow;

ax_pos = ax_main.Position;
xlim_main = ax_main.XLim;
ylim_main = ax_main.YLim;

% Centro del lato superiore del rettangolo
x_center = ax_pos(1) + ...
    ((1 + 7)/2 - xlim_main(1)) / diff(xlim_main) * ax_pos(3);

y_top = ax_pos(2) + ...
    ((ymax + y_margin) - ylim_main(1)) / diff(ylim_main) * ax_pos(4);

% Centro del lato inferiore dell'inset
inset_pos = ax_inset.Position;
x_inset = inset_pos(1) + inset_pos(3)/2;
y_inset = inset_pos(2);

annotation('line', ...
    [x_center x_inset], ...
    [y_top y_inset], ...
    'LineStyle', '-', ...
    'Color', 'k', ...
    'LineWidth', 0.5);