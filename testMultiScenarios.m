clear; close all; clc;

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
N0 = 1e5;
I0 = 70;
A0 = 30;


t = table;
t.dE2A = linspace(3, 10, 8)';
t.dE2I = linspace(3, 10, 8)';
t.dI2R = linspace(3, 10, 8)';
t.dI2H = linspace(3, 10, 8)';
t.pE2A = linspace(0.1, 0.8, 8)';
t.pI2R = linspace(0.1, 0.8, 8)';

params0 = params;

indicesTable = table;
for k = 1:size(t,2)
    for i = 1:size(t,1)
        S = struct;
        S.type = '.';
        S.subs = t.Properties.VariableNames{k};
        values = subsref(t, S);
        value_i = values(i);
        params = subsasgn(params0, S, value_i);

        md1 = OrigionalModel('Params', params, 'x0', [N0, 0, I0, A0, 0, 0]');
        md2 = ImprovedModel('Params', params, 'x0', [N0, 0, 0, I0 * params.pI2R, I0 * (1-params.pI2R), A0, 0, 0]');
        md3 = ComplexModel('Params', params, "N", N0, "I", I0, 'A', A0);

        solve(md1);
        solve(md2);

        tic;
        solve(md3);
        toc

        %% compute and save indices
        tempTable = table;
        tempTable.dE2A = params.dE2A;
        tempTable.dE2I = params.dE2I;
        tempTable.dI2R = params.dI2R;
        tempTable.dI2H = params.dI2H;
        tempTable.pE2A = params.pE2A;
        tempTable.pI2R = params.pI2R;

        % A(t)
        [tempTable.maxA1, idxA1] = max(md1.record.A);
        [tempTable.maxA2, idxA2] = max(md2.record.A);
        [tempTable.maxA3, idxA3] = max(md3.record.A);
        tempTable.tA1 = md1.record.tGrids(idxA1);
        tempTable.tA2 = md2.record.tGrids(idxA2);
        tempTable.tA3 = md3.record.tGrids(idxA3);
        
        % H(t)
        [tempTable.maxH1, idxH1] = max(md1.record.H);
        [tempTable.maxH2, idxH2] = max(md2.record.H);
        [tempTable.maxH3, idxH3] = max(md3.record.H);
        tempTable.tH1 = md1.record.tGrids(idxH1);
        tempTable.tH2 = md2.record.tGrids(idxH2);
        tempTable.tH3 = md3.record.tGrids(idxH3);

        % infectionRate(t)
        [tempTable.maxInfectionRates1, idxInfectionRates1] = max(md1.record.InfectionRates);
        [tempTable.maxInfectionRates2, idxInfectionRates2] = max(md2.record.InfectionRates);
        [tempTable.maxInfectionRates3, idxInfectionRates3] = max(md3.record.InfectionRates);
        tempTable.tInfectionRates1 = md1.record.tGrids(idxInfectionRates1);
        tempTable.tInfectionRates2 = md2.record.tGrids(idxInfectionRates2);
        tempTable.tInfectionRates3 = md3.record.tGrids(idxInfectionRates3);

        % cumulative infection trapz(infectionRate)
        tempTable.cumulativeInfection1 = trapz(md1.record.tGrids, md1.record.InfectionRates);
        tempTable.cumulativeInfection2 = trapz(md2.record.tGrids, md2.record.InfectionRates);
        tempTable.cumulativeInfection3 = sum(md3.record.InfectionRates);
        
        indicesTable = [indicesTable; tempTable];

        %%
        fig = figure;
        fig.Position = [1891,-288,1116,1283];
        tilde = tiledlayout(4,2);

        nexttile;
        plot(md1.record.tGrids, md1.record.S, 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.S, 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.S, 'g-', 'LineWidth', 2);
        %legend('Classic Model', 'Improved Model', 'Complex Model');
        ylabel('$$S(t)$$', 'Interpreter', 'latex')

        nexttile;
        plot(md1.record.tGrids, md1.record.E, 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.E, 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.E, 'g-', 'LineWidth', 2);
        %legend('Classic Model', 'Improved Model', 'Complex Model');
        ylabel('$$E(t)$$', 'Interpreter', 'latex')

        nexttile;
        plot(md1.record.tGrids, md1.record.I, 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.I, 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.I, 'g-', 'LineWidth', 2);
        %legend('Classic Model', 'Improved Model', 'Complex Model');
        ylabel('$$I(t)$$', 'Interpreter', 'latex')

        nexttile;
        plot(md1.record.tGrids, md1.record.R, 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.R, 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.R, 'g-', 'LineWidth', 2);
        %legend('Classic Model', 'Improved Model', 'Complex Model');
        ylabel('$$R(t)$$', 'Interpreter', 'latex')


        nexttile;
        plot(md1.record.tGrids, md1.record.A, 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.A, 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.A, 'g-', 'LineWidth', 2);
        %legend('Classic Model', 'Improved Model', 'Complex Model');
        ylabel('$$A(t)$$', 'Interpreter', 'latex')

        nexttile;
        plot(md1.record.tGrids, md1.record.H, 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.H, 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.H, 'g-', 'LineWidth', 2);
        %legend('Classic Model', 'Improved Model', 'Complex Model');
        ylabel('$$H(t)$$', 'Interpreter', 'latex')

        ax = nexttile;
        plot(md1.record.tGrids, md1.record.CumulativeA ./ (md1.record.CumulativeA + md1.record.CumulativeI), 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, md2.record.CumulativeA ./ (md2.record.CumulativeA + md2.record.CumulativeI), 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, md3.record.CumulativeA ./ (md3.record.CumulativeA + md3.record.CumulativeI), 'g-', 'LineWidth', 2);
        plot(md1.record.tGrids, params.pE2A * ones(size(md1.record.tGrids)), 'black-.', 'LineWidth', 2);
        %ylabel('$$\sum A / (\sum A + \sum I)$$', 'Interpreter', 'latex');
        ylabel('Asymptomatic Ratio');
        %legend('Classic Model', 'Improved Model', 'Complex Model', 'Real Value', 'Location', 'southeast');
        ax.YLim(1) = ax.YLim(1) * 0.8;
        ax.YLim(2) = ax.YLim(2) * 1.2;

        nexttile;
        plot(md1.record.tGrids, cumtrapz(md1.params.pI2R*md1.record.I/md1.params.dI2R) ./ cumtrapz(md1.params.pI2R*md1.record.I/md1.params.dI2R + md1.record.H./md1.params.dH2R), 'b-', 'LineWidth', 2); hold on;
        plot(md2.record.tGrids, cumtrapz(md2.record.I1/md2.params.dI2R) ./ cumtrapz(md2.record.I1/md2.params.dI2R + md2.record.H./md2.params.dH2R), 'r-', 'LineWidth', 2);
        plot(md3.record.tGrids, cumsum(md3.record.newI2R) ./ cumsum(md3.record.newI2R + md3.record.newH2R), 'g-', 'LineWidth', 2);
        plot(md1.record.tGrids, params.pI2R * ones(size(md1.record.tGrids)), 'black-.', 'LineWidth', 2);
        %ylabel('$$\sum R / (\sum R + \sum C)$$', 'Interpreter', 'latex');
        ylabel('Self-recovered Ratio');
        %legend('Classic Model', 'Improved Model', 'Complex Model', 'Real Value');


        xlabel(tilde, 'Time (in days)', 'FontName', 'times new roman', 'FontSize', 18);

        for ii = 1:numel(tilde.Children)
            if isa(tilde.Children(ii), "matlab.graphics.axis.Axes")
                tilde.Children(ii).FontName = 'times new roman';
                tilde.Children(ii).FontSize = 18;
            end
        end


        % 创建 textbox
        annotation(fig,'textbox',[0.0708924731182796 0.943881527669525 0.059931899641577 0.029618082618862],'String','A.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.0735806451612903 0.271239282930632 0.059931899641577 0.0296180826188621],'String','G.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.511752688172043 0.273577552611068 0.0599318996415771 0.0296180826188621],'String','H.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.498311827956989 0.9423226812159 0.0599318996415771 0.029618082618862],'String','B.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.0753727598566307 0.709275136399065 0.0599318996415769 0.029618082618862],'String','C.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.514440860215054 0.708495713172253 0.0599318996415771 0.029618082618862],'String','D.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.0708924731182796 0.489477786438036 0.059931899641577 0.0296180826188621],'String','E.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % 创建 textbox
        annotation(fig,'textbox',[0.510856630824373 0.493374902572097 0.0599318996415771 0.029618082618862],'String','F.','FontSize',18,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');


        exportgraphics(fig, "Result_" + S.subs + "=" + value_i + ".jpg", 'Resolution', 300);


        close all;
    end
end

writetable(indicesTable, 'indicesTable.xlsx');