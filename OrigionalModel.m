% SEIARH
classdef OrigionalModel < handle

    properties
        % model parameters
        params

        % initial values
        x0 {mustBeNonnegative}

        % 
        record
    end

    methods

        function md = OrigionalModel(args)
            arguments
                args.Params = []; 
                args.x0 = [];
            end

            if ~isempty(args.Params)
                md.params = args.Params;
            end

            if ~isempty(args.x0)
                md.x0 = args.x0;
            end

        end

        function [dxdt, newI, newA, newR, newH] = computeDerivative(md, t, x)
            % x = [S,E,I,A,R,H]
            S = x(1);
            E = x(2);
            I = x(3);
            A = x(4);
            R = x(5);
            H = x(6);

            rateS2E = md.params.beta * S * (I + md.params.kappa * A) / sum(x);
            rateE2A = md.params.pE2A * E / md.params.dE2A;
            rateE2I = (1-md.params.pE2A) * E / md.params.dE2I;

            rateA2R = A / md.params.dA2R;
            rateI2R = md.params.pI2R * I / md.params.dI2R;
            rateI2H = (1-md.params.pI2R) * I / md.params.dI2H;
            rateH2R = H / md.params.dH2R;

            dSdt = -rateS2E;
            dEdt = rateS2E - rateE2A - rateE2I;
            dIdt = rateE2I - rateI2R - rateI2H;
            dAdt = rateE2A - rateA2R;
            dRdt = rateA2R + rateI2R + rateH2R;
            dHdt = rateI2H - rateH2R;

            newI = rateE2I;
            newA = rateE2A;
            newR = rateI2R +  rateA2R + rateH2R;
            newH = rateI2H;
            dxdt = [dSdt; dEdt; dIdt; dAdt; dRdt; dHdt];
        end

        function solve(md)
            [tt, xx] = ode45(@(t,x) md.computeDerivative(t,x), md.params.tGrids, md.x0);
            variableNames = {'S','E','I','A','R','H'};

            temp = table;
            temp.tGrids = tt;

            md.record = [temp, array2table(xx, 'VariableNames', variableNames)];
            
            % save infection rates
            dxdt = zeros(size(xx));
            newI = zeros(numel(tt), 1);
            newA = zeros(numel(tt), 1);
            newR = zeros(numel(tt), 1);
            newH = zeros(numel(tt), 1);
            for i = 1:numel(tt)
                [dxdt(i,:), newI(i), newA(i), newR(i), newH(i)] = md.computeDerivative(tt(i), xx(i,:)');
            end
            md.record.InfectionRates = -dxdt(:,1);
            md.record.CumulativeI = cumtrapz(tt, newI);
            md.record.CumulativeA = cumtrapz(tt, newA);
            md.record.CumulativeR = cumtrapz(tt, newR);
            md.record.CumulativeH = cumtrapz(tt, newH);
        end
    end
end