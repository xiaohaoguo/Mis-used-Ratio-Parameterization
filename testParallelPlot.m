clear; close all; clc;

dataTable = readtable("indicesTable.xlsx");
dataTable(:,1:6) = [];

m = size(dataTable, 1); % number of rows
n = size(dataTable, 2) / 3; % number of indices
idxn = 1:3:(3*n-2);
yData = [dataTable{:,idxn}; dataTable{:,idxn+1}; dataTable{:,idxn+2}];

%GroupBy = kron((1:3)', ones(m,1));
idxm = 1:3:(3*m-2);
xGroupData = repmat(["Classic Model"; "Improved Model"; "Complex Model"], m, 1);
xGroupData = categorical(xGroupData([idxm, idxm+1, idxm+2]));

variableNames1 = {'max A(t)', 'argmax A(t)', 'max H(t)', 'argmax H(t)', 'max InfectionRate(t)', 'argmax InfectionRate(t)', 'Cumulative Infection'};
variableNames2 = {'max $A(t)$', 'argmax $A(t)$', 'max $H(t)$', 'argmax $H(t)$', 'max $InfectionRate(t)$', 'argmax $InfectionRate(t)$', 'Cumulative Infection'};

fig = figure;
fig.WindowState = 'maximized';
tilde = tiledlayout(2,7,'TileSpacing','tight');

ax = nexttile([1,7]);
p = parallelplot(yData,'GroupData',xGroupData, 'DataNormalization', 'norm');
%legend({'Classic Model', 'Improved Model', 'Complex Model'});
p.CoordinateTickLabels = variableNames1;
set(gca, 'FontSize', 20, 'FontName', 'times new roman');
%p.CoordinateLabel

% boxchart for 7 indices
for i = 1:7 
    ax = nexttile;
    b = boxchart(xGroupData, yData(:,i), 'Notch', 'on');
    %title(gca, variableNames2(i), 'Interpreter', 'latex');
    set(gca, 'FontSize', 16, 'FontName', 'times new roman');
end



annotation(fig,'textbox',...
    [0.0263906250000001 0.920091324200911 0.020484375 0.0372907153729072],...
    'String','A.',...
    'FontSize',22,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(fig,'textbox',...
    [0.0275625000000001 0.489345509893454 0.020484375 0.0372907153729072],...
    'String','B.',...
    'FontSize',22,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');


exportgraphics(fig, 'Boxchart.jpg', 'resolution', 300);