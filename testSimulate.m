clear; close all; clc;

params = struct;
params.beta = 0.3;
params.c = 10;
params.q = params.beta / params.c;
params.kappa = 1;
params.pE2A = 0.3;
params.pI2R = 0.8;
%params.pI2R = 0.2;

params.dE2A = 5;
params.dE2I = 3;

params.dI2R = 7;
params.dI2H = 3;
params.dA2R = 5;
params.dH2R = 11;

params.tGrids = linspace(0, 300, 1e3);
N0 = 1e5;
I0 = 70;
A0 = 30;

md1 = OrigionalModel('Params', params, 'x0', [N0, 0, I0, A0, 0, 0]');
md2 = ImprovedModel('Params', params, 'x0', [N0, 0, 0, I0 * params.pI2R, I0 * (1-params.pI2R), A0, 0, 0]');
md3 = ComplexModel('Params', params, "N", N0, "I", I0, 'A', A0);

solve(md1);
solve(md2);

tic; 
solve(md3); 
toc

%%
fig = figure; 
fig.Position = [407,230,1193,843];
tilde = tiledlayout(2,2);

nexttile;
plot(md1.record.tGrids, md1.record.A, 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.A, 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.A, 'g-', 'LineWidth', 2);
legend('Classic Model', 'Improved Model', 'Complex Model');
ylabel('$$A(t)$$', 'Interpreter', 'latex')

nexttile;
plot(md1.record.tGrids, md1.record.H, 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.H, 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.H, 'g-', 'LineWidth', 2);
legend('Classic Model', 'Improved Model', 'Complex Model');
ylabel('$$H(t)$$', 'Interpreter', 'latex')

ax = nexttile;
plot(md1.record.tGrids, md1.record.CumulativeA ./ (md1.record.CumulativeA + md1.record.CumulativeI), 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.CumulativeA ./ (md2.record.CumulativeA + md2.record.CumulativeI), 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.CumulativeA ./ (md3.record.CumulativeA + md3.record.CumulativeI), 'g-', 'LineWidth', 2);
plot(md1.record.tGrids, params.pE2A * ones(size(md1.record.tGrids)), 'black-.', 'LineWidth', 2);
%ylabel('$$\sum A / (\sum A + \sum I)$$', 'Interpreter', 'latex');
ylabel('Asymptomatic Ratio');
legend('Classic Model', 'Improved Model', 'Complex Model', 'Real Value', 'Location', 'southeast');
ax.YLim(1) = ax.YLim(1) * 0.8;
ax.YLim(2) = ax.YLim(2) * 1.2;

nexttile;
plot(md1.record.tGrids, cumtrapz(md1.params.pI2R*md1.record.I/md1.params.dI2R) ./ cumtrapz(md1.params.pI2R*md1.record.I/md1.params.dI2R + md1.record.H./md1.params.dH2R), 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, cumtrapz(md2.record.I1/md2.params.dI2R) ./ cumtrapz(md2.record.I1/md2.params.dI2R + md2.record.H./md2.params.dH2R), 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, cumsum(md3.record.newI2R) ./ cumsum(md3.record.newI2R + md3.record.newH2R), 'g-', 'LineWidth', 2);
plot(md1.record.tGrids, params.pI2R * ones(size(md1.record.tGrids)), 'black-.', 'LineWidth', 2);
%ylabel('$$\sum R / (\sum R + \sum C)$$', 'Interpreter', 'latex');
ylabel('Self-recovered Ratio');
legend('Classic Model', 'Improved Model', 'Complex Model', 'Real Value');


xlabel(tilde, 'Time (in days)', 'FontName', 'times new roman', 'FontSize', 18);

for i = 1:numel(tilde.Children)
    if isa(tilde.Children(i), "matlab.graphics.axis.Axes")
        tilde.Children(i).FontName = 'times new roman';
        tilde.Children(i).FontSize = 18;
    end
end

annotation(fig,'textbox',[0.0316946688206787 0.931318681318683 0.0345411954765751 0.0631868131868131],'String',{'A.'},'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% 创建 textbox
annotation(fig,'textbox',[0.492114701130857 0.928571428571431 0.034541195476575 0.0631868131868131],'String','B.','FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% 创建 textbox
annotation(fig,'textbox',[0.0325024232633282 0.457417582417586 0.0345411954765751 0.0631868131868131],'String','C.','FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% 创建 textbox
annotation(fig,'textbox',[0.495345718901456 0.454670329670333 0.034541195476575 0.0631868131868131],'String','D.','FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

exportgraphics(fig, 'NumericalExperiment.jpg', 'Resolution', 300);



fig2 = figure;
fig2.Position = [1308,145,1133,855];
tilde2 = tiledlayout(2,2);


nexttile;
plot(md1.record.tGrids, md1.record.S, 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.S, 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.S, 'g-', 'LineWidth', 2);
legend('Classic Model', 'Improved Model', 'Complex Model');
ylabel('$$S(t)$$', 'Interpreter', 'latex')

nexttile;
plot(md1.record.tGrids, md1.record.E, 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.E, 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.E, 'g-', 'LineWidth', 2);
legend('Classic Model', 'Improved Model', 'Complex Model');
ylabel('$$E(t)$$', 'Interpreter', 'latex')

nexttile;
plot(md1.record.tGrids, md1.record.I, 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.I, 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.I, 'g-', 'LineWidth', 2);
legend('Classic Model', 'Improved Model', 'Complex Model');
ylabel('$$I(t)$$', 'Interpreter', 'latex')

nexttile;
plot(md1.record.tGrids, md1.record.R, 'b-', 'LineWidth', 2); hold on;
plot(md2.record.tGrids, md2.record.R, 'r-', 'LineWidth', 2);
plot(md3.record.tGrids, md3.record.R, 'g-', 'LineWidth', 2);
legend('Classic Model', 'Improved Model', 'Complex Model');
ylabel('$$R(t)$$', 'Interpreter', 'latex')

for i = 1:numel(tilde2.Children)
    if isa(tilde2.Children(i), "matlab.graphics.axis.Axes")
        tilde2.Children(i).FontName = 'times new roman';
        tilde2.Children(i).FontSize = 18;
    end
end