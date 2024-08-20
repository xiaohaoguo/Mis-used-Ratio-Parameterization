% SEIARH
classdef ImprovedModel < handle

    properties
        % model parameters
        params

        % initial values
        x0 {mustBeNonnegative}

        % 
        record
    end

    methods

        function md = ImprovedModel(args)
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
            % x = [S,E1,E2,I1,I2,A,R,H]
            S = x(1);
            E1 = x(2);
            E2 = x(3);
            I1 = x(4);
            I2 = x(5);
            A = x(6);
            R = x(7);
            H = x(8);

            I = I1 + I2;

            rateS2E1 = md.params.pE2A * md.params.beta * S * (I + md.params.kappa * A) / sum(x);
            rateS2E2 = (1-md.params.pE2A) * md.params.beta * S * (I + md.params.kappa * A) / sum(x);
            rateE2A = E1 / md.params.dE2A;
            rateE2I1 = md.params.pI2R * E2 / md.params.dE2I;
            rateE2I2 = (1-md.params.pI2R) * E2 / md.params.dE2I;

            rateA2R = A / md.params.dA2R;
            rateI2R = I1 / md.params.dI2R;
            rateI2H = I2 / md.params.dI2H;
            rateH2R = H / md.params.dH2R;

            dSdt = -rateS2E1 - rateS2E2;
            dE1dt = rateS2E1 - rateE2A;
            dE2dt = rateS2E2 - rateE2I1 - rateE2I2;
            dI1dt = rateE2I1 - rateI2R;
            dI2dt = rateE2I2 - rateI2H;
            dAdt = rateE2A - rateA2R;
            dRdt = rateA2R + rateI2R + rateH2R;
            dHdt = rateI2H - rateH2R;

            newI = rateE2I1 + rateE2I2;
            newA = rateE2A;
            newR = rateI2R + rateA2R + rateH2R;
            newH = rateI2H;
            dxdt = [dSdt; dE1dt; dE2dt; dI1dt; dI2dt; dAdt; dRdt; dHdt];
        end

        function solve(md)
            [tt, xx] = ode45(@(t,x) md.computeDerivative(t,x), md.params.tGrids, md.x0);
            variableNames = {'S','E1', 'E2','I1','I2','A','R','H'};

            temp = table;
            temp.tGrids = tt;

            md.record = [temp, array2table(xx, 'VariableNames', variableNames)];
            md.record.E = md.record.E1 + md.record.E2;
            md.record.I = md.record.I1 + md.record.I2;
            md.record.I1 = md.record.I1;

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