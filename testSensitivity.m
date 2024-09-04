clear; close all; clc;

% params.pE2A = 0.3;
% params.pI2R = 0.8;
N0 = 1e5;
I0 = 70;
A0 = 30;

dI2R = 1 : 0.1 : 10;
dI2H = 1 : 0.1 : 10;
dE2A = 1 : 0.1 : 10;
dE2I = 1 : 0.1 : 10;

m = numel(dI2R);
n = numel(dI2H);


cumulativeH1 = zeros(numel(dI2R), numel(dI2H));
cumulativeH2 = zeros(numel(dI2R), numel(dI2H));

cumulativeA1 = zeros(numel(dE2A), numel(dE2I));
cumulativeA2 = zeros(numel(dE2A), numel(dE2I));

%%% for over/under estimate of Hospitalizations
tic;
parfor i = 1:numel(dI2R)
    for j = 1:n %numel(dI2H)
        params = initializeParameters();
        params.dI2R = dI2R(i);
        params.dI2H = dI2H(j);

        md1 = OrigionalModel('Params', params, 'x0', [N0, 0, I0, A0, 0, 0]');
        md2 = ImprovedModel('Params', params, 'x0', [N0, 0, 0, I0 * params.pI2R, I0 * (1-params.pI2R), A0, 0, 0]');
        solve(md1);
        solve(md2);

        % save the 'Cumulative H' simulated by both models
        cumulativeH1(i,j) = trapz(md1.record.tGrids, (1-params.pI2R)/params.dI2H*md1.record.I);
        cumulativeH2(i,j) = trapz(md2.record.tGrids, md2.record.I2/params.dI2H);
    end
end
toc;


%%% for over/under estimate of Asymptomatics
tic;
parfor i = 1:numel(dE2A)
    for j = 1:n %numel(dE2I)
        params = initializeParameters();
        params.dE2A = dE2A(i);
        params.dE2I = dE2I(j);

        md1 = OrigionalModel('Params', params, 'x0', [N0, 0, I0, 0, 0, 0]');
        md2 = ImprovedModel('Params', params, 'x0', [N0, 0, 0, I0 * params.pI2R, I0 * (1-params.pI2R), 0, 0, 0]');
        solve(md1);
        solve(md2);

        % save the 'Cumulative H' simulated by both models
        cumulativeA1(i,j) = trapz(md1.record.tGrids, params.pE2A/params.dE2A*md1.record.E);
        cumulativeA2(i,j) = trapz(md2.record.tGrids, md2.record.E1/params.dE2A);
    end
end
toc;


%%
ampH = (cumulativeH1 - cumulativeH2) ./ cumulativeH2;
ampA = (cumulativeA1 - cumulativeA2) ./ cumulativeA2;

[X, Y] = meshgrid(dI2R, dI2H);
XX = linspace(min(dI2R), max(dI2R), 1e3);
YY = linspace(min(dI2H), max(dI2H), 1e3);
[XX, YY] = meshgrid(XX,YY);

ampH = interp2(X, Y, ampH, XX, YY);
ampA = interp2(X, Y, ampA, XX, YY);

ampUpperH = ampH;
ampUpperH(ampH < 0) = nan;
ampLowerH = ampH;
ampLowerH(ampH > 0) = nan;

ampUpperA = ampA;
ampUpperA(ampA < 0) = nan;
ampLowerA = ampA;
ampLowerA(ampA > 0) = nan;


fig3 = figure;
fig3.Position = [327,578,1233,660];
tiled3 = tiledlayout(1,2);

ax5 = nexttile;
levels1 = [-inf, -1 : 0.2 : -0.2, 0.2 : 0.2 : 1.4, inf];
contourf(XX, YY, ampH, levels1, 'ShowText','on'); hold on;
contour(XX, YY, ampH, [0,0], 'ShowText','on', 'LineWidth', 2, 'LineStyle', '-', 'EdgeColor', 'black');

cbar = colorbar;
colormap('sky');
clim([-1, 1.4]);
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'times new roman');
set(gca, 'XTick', 1:10);
set(gca, 'YTick', 1:10);
axis equal;
title(ax5, 'Relative bias for hospitalization', 'FontSize', 18, 'FontName', 'times new roman')

ax6 = nexttile;
levels2 = [-0.8 : 0.2 :-0.2, 0.2 : 0.2 : 1.2];
contourf(XX, YY, ampA, levels2, 'ShowText','on'); hold on;
contour(XX, YY, ampA, [0,0], 'ShowText','on', 'LineWidth', 2, 'LineStyle', '-', 'EdgeColor', 'black');
colorbar;
colormap('sky')
clim([-0.8, 1.2]);
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'times new roman');
set(gca, 'XTick', 1:10);
set(gca, 'YTick', 1:10);
axis equal;
title(ax6, 'Relative bias for asymptomatics', 'FontSize', 18, 'FontName', 'times new roman')

xlabel(ax5, "Period from I to R ($1/\gamma$, in days)", 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'times new roman');
ylabel(ax5, "Period from I to H ($1/\gamma^{'}$, in days)", 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'times new roman');
xlabel(ax6, "Period from E to A ($1/\omega$, in days)", 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'times new roman');
ylabel(ax6, "Period from E to I ($1/\omega^{'}$, in days)", 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'times new roman');


% 创建 textbox
annotation(fig3,'textbox',...
    [0.0804809407948094 0.837878787878788 0.0922684509326846 0.093939393939394],...
    'String',{'A.'},...
    'FontSize',24,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% 创建 textbox
annotation(fig3,'textbox',...
    [0.546012165450122 0.848484848484849 0.0922684509326848 0.0939393939393939],...
    'String','B.',...
    'FontSize',24,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

exportgraphics(fig3, "RelativeBias.jpg", 'Resolution', 300);

function params = initializeParameters()
    params = struct;
    params.beta = 0.3;
    params.c = 10;
    params.q = params.beta / params.c;
    params.kappa = 1;
    params.pE2A = 0.3;
    params.pI2R = 0.8;
    
    params.dE2A = 5;
    params.dE2I = 3;
    
    params.dI2R = 7;
    params.dI2H = 3;
    
    params.dA2R = 5;
    params.dH2R = 11;
    
    params.tGrids = linspace(0, 300, 1e3);
end