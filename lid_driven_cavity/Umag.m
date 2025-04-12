clear
close all
% format compact
% set(groot,'defaulttextinterpreter','latex');
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');

data = load('Umag.txt');
x = data(:,1);
y = data(:,2);
V = data(:,3);

% Definir grilla regular
xlin = linspace(min(x), max(x), 100);  % puedes ajustar la resolución (100, 200, etc.)
ylin = linspace(min(y), max(y), 100);
[X, Y] = meshgrid(xlin, ylin);

% Interpolar los valores de velocidad sobre la grilla
Vgrid = griddata(x, y, V, X, Y, 'linear');  % puedes usar 'linear', 'nearest', 'cubic'

contourf(X, Y, Vgrid, 20);  % 20 líneas de nivel, ajustable
xlabel('$x$','FontSize', 14, 'Interpreter', 'latex')
ylabel('$y$','FontSize', 14, 'Interpreter', 'latex')
title('Magnitud de la velocidad','FontSize', 16, 'Interpreter', 'latex')
set(gca, 'FontSize', 12);  % tamaño de los números en los ejes
colorbar
caxis([0 1]);
colormap jet
axis equal
