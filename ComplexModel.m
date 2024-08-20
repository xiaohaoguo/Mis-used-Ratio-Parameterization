% SEIARH
classdef ComplexModel < handle

    properties
        %%%Individual Properties
        % Disease stage of individuals
        % 1:6 represent SEIARH respectively
        % 
        % default: ones(N,1)
        stage (:,1)

        % Time elapsed since entered the current stage
        duration (:,1)

        % random setting of individual transitions
        idxE2A
        idxI2R

        %%% parameters
        % StepSize for simulation
        dt = 0.1

        % struct of other parameters
        params

        %%% simulated records
        record
    end

    methods
        
        function md = ComplexModel(args)
            arguments
                args.Params = [];
                args.N = 1e5
                args.I = 70
                args.A = 30
            end

            if ~isempty(args.Params)
                md.params = args.Params;
            end

            md.stage = ones(args.N, 1);
            md.stage(1:args.I) = 3;
            md.stage(args.I+1 : args.I +args.A) = 4;
            md.duration = zeros(args.N,1);
        end

        function info = updateStages(md)
            % SEIARH

            idxS = md.stage == 1;
            idxE = md.stage == 2;
            idxI = md.stage == 3;
            idxA = md.stage == 4;
            idxR = md.stage == 5;
            idxH = md.stage == 6;

            %% transition
            tI = md.dt * ones(size(md.stage));
            tA = md.dt * ones(size(md.stage));


            % H to R
            idxDevelopedH2R = md.duration + md.dt >= md.params.dH2R;
            idxH2Red = idxH & idxDevelopedH2R;
            idxH2Ring = idxH & ~idxDevelopedH2R;
            md.stage(idxH2Red) = 5;
            md.duration(idxH2Red) = md.duration(idxH2Red) + md.dt - md.params.dH2R(idxH2Red);
            md.duration(idxH2Ring) = md.duration(idxH2Ring) + md.dt;
            
            % I to H/R
            idxDevelopedI2R = md.duration + md.dt >= md.params.dI2R;
            idxI2Red = idxI & md.idxI2R & idxDevelopedI2R;
            idxI2Ring= idxI & md.idxI2R & ~idxDevelopedI2R;
            md.stage(idxI2Red) = 5;
            md.duration(idxI2Red) = md.duration(idxI2Red) + md.dt - md.params.dI2R(idxI2Red);
            md.duration(idxI2Ring) = md.duration(idxI2Ring) + md.dt;
            tI(idxI2Red) = md.dt - md.duration(idxI2Red); % half-way transmissions



            idxDevelopedI2H = (md.duration + md.dt >= md.params.dI2H);
            idxI2Hed = idxI & (~md.idxI2R & ~md.idxE2A) & idxDevelopedI2H & (md.duration + md.dt < md.params.dI2H + md.params.dH2R);
            idxI2Hing= idxI & (~md.idxI2R & ~md.idxE2A) & ~idxDevelopedI2H;
            idxI2H2Red = idxI & (~md.idxI2R & ~md.idxE2A) & (md.duration + md.dt >= md.params.dI2H + md.params.dH2R);
            md.stage(idxI2Hed) = 6;
            md.stage(idxI2H2Red) = 5;
            md.duration(idxI2Hed) = md.duration(idxI2Hed) + md.dt - md.params.dI2H(idxI2Hed);
            md.duration(idxI2Hing) = md.duration(idxI2Hing) + md.dt;
            md.duration(idxI2H2Red) = md.duration(idxI2H2Red) + md.dt - (md.params.dI2H(idxI2H2Red) + md.params.dH2R(idxI2H2Red));
            tI(idxI2Hed) = md.dt - md.duration(idxI2Hed); % half-way transmissions
            tI(idxI2H2Red) = md.params.dI2H(idxI2H2Red);

 
            % A to R
            idxDevelopedA2R = md.duration + md.dt >= md.params.dA2R;
            idxA2Red = idxA & idxDevelopedA2R;
            idxA2Ring = idxA & ~idxDevelopedA2R;
            md.stage(idxA2Red) = 5;
            md.duration(idxA2Red) = md.duration(idxA2Red) + md.dt - md.params.dA2R(idxA2Red);
            md.duration(idxA2Ring) = md.duration(idxA2Ring) + md.dt;
            tA(idxA2Red) = md.dt - md.duration(idxA2Red); % half-way transmissions
             

            % E to A/I
            idxDevelopedE2A = md.duration + md.dt >= md.params.dE2A;
            idxE2Aed = idxE & md.idxE2A & idxDevelopedE2A & (md.duration + md.dt < md.params.dE2A + md.params.dA2R);
            idxE2Aing= idxE & md.idxE2A & ~idxDevelopedE2A;
            idxE2A2Red = idxE & md.idxE2A & (md.duration + md.dt >= md.params.dE2A + md.params.dA2R);
            md.stage(idxE2Aed) = 4;
            md.stage(idxE2A2Red) = 5;
            md.duration(idxE2Aed) = md.duration(idxE2Aed) + md.dt - md.params.dE2A(idxE2Aed);
            md.duration(idxE2Aing) = md.duration(idxE2Aing) + md.dt;
            md.duration(idxE2A2Red) = md.duration(idxE2A2Red) + md.dt - (md.params.dE2A(idxE2A2Red) + md.params.dA2R(idxE2A2Red));
            tA(idxE2Aed) = md.duration(idxE2Aed); % half-way transmissions
            tA(idxE2A2Red) = md.params.dA2R(idxE2A2Red);
            
            

            idxDevelopedE2I = md.duration + md.dt >= md.params.dE2I;
            idxE2Ied = idxE & ~md.idxE2A & idxDevelopedE2I & (  ((md.duration + md.dt < md.params.dE2I + md.params.dI2H) & ~md.idxI2R)...
                                                                         | ((md.duration + md.dt < md.params.dE2I + md.params.dI2R) & md.idxI2R) );
            idxE2Iing = idxE & ~md.idxE2A & ~idxDevelopedE2I;
            idxE2I2Hed = idxE & ~md.idxE2A & ~md.idxI2R & (md.duration + md.dt >= md.params.dE2I + md.params.dI2H)...
                                                        & (md.duration + md.dt < md.params.dE2I + md.params.dI2H + md.params.dH2R);
            idxE2I2Red = idxE & ~md.idxE2A & md.idxI2R & (md.duration + md.dt >= md.params.dE2I + md.params.dI2R);
            idxE2I2H2Red = idxE & ~md.idxE2A & ~md.idxI2R & (md.duration + md.dt >= md.params.dE2I + md.params.dI2H + md.params.dH2R);
            md.stage(idxE2Ied) = 3;
            md.stage(idxE2I2Hed) = 6;
            md.stage(idxE2I2Red) = 5;
            md.stage(idxE2I2H2Red) = 5;
            md.duration(idxE2Ied) = md.duration(idxE2Ied) + md.dt - md.params.dE2I(idxE2Ied);
            md.duration(idxE2Iing) = md.duration(idxE2Iing) + md.dt;
            md.duration(idxE2I2Hed) = md.duration(idxE2I2Hed) + md.dt - (md.params.dE2I(idxE2I2Hed) + md.params.dI2H(idxE2I2Hed));
            md.duration(idxE2I2Red) = md.duration(idxE2I2Red) + md.dt - (md.params.dE2I(idxE2I2Red) + md.params.dI2R(idxE2I2Red));
            md.duration(idxE2I2H2Red) = md.duration(idxE2I2H2Red) + md.dt - (md.params.dE2I(idxE2I2H2Red) + md.params.dI2H(idxE2I2H2Red) + md.params.dH2R(idxE2I2H2Red));
            tI(idxE2Ied) = md.duration(idxE2Ied); % half-way transmissions
            tI(idxE2I2Hed | idxE2I2H2Red) = md.params.dI2R(idxE2I2Hed | idxE2I2H2Red);
            tI(idxE2I2Red) = md.params.dI2R(idxE2I2Red);


            if any(md.duration < 0)
                1;
            end

            % for the case that dE2I + dI2R < dt, that tI < 0
            % tI = max(0, tI);
            % tA = max(0, tA);

            %% Transmission
            cMat = randi(numel(md.stage), [numel(md.stage), md.params.c]);

            % probability of infection (for each individuals)
            tempI = tI .* idxI;
            tempA = tA .* idxA;
            Prob = md.params.q .* sum(tempI(cMat) + md.params.kappa*tempA(cMat), 2); 
            %Prob = md.params.q .* sum(idxI(cMat) + md.params.kappa*idxA(cMat), 2); 
            %Prob = 1 - (1 - md.params.q) .^ sum(tempI(cMat) + md.params.kappa*tempA(cMat), 2); 
            %Prob = 1 - (1 - md.params.q) .^ (md.dt * sum(idxI(cMat) + md.params.kappa*idxA(cMat), 2)); 
            
            isInfected = idxS & (rand(numel(md.stage),1) < Prob);
            md.stage(isInfected) = 2;
            md.duration(isInfected) = 0;

            info = table;
            info.S = sum(md.stage == 1);
            info.E = sum(md.stage == 2);
            info.I = sum(md.stage == 3);
            info.A = sum(md.stage == 4);
            info.R = sum(md.stage == 5);
            info.H = sum(md.stage == 6);
            info.InfectionRates = sum(isInfected);
            info.newI = sum(idxE2Ied);
            info.newA = sum(idxE2Aed);
            info.newR = sum(idxI2Red | idxA2Red | idxH2Red);
            info.newH = sum(idxI2Hed);
            info.newI2R = sum(idxI2Red);
            info.newH2R = sum(idxH2Red);

        end %updateStage

        function solve(md)
            % initialize random settings of individual transition
            md.idxE2A = rand(size(md.stage)) < md.params.pE2A;
            md.idxI2R = rand(size(md.stage)) < md.params.pI2R;

            % randomize params.d*
            md.params.dE2I = replaceConstantByDistribution(md.params.dE2I, numel(md.stage));
            md.params.dE2A = replaceConstantByDistribution(md.params.dE2A, numel(md.stage));
            md.params.dA2R = replaceConstantByDistribution(md.params.dA2R, numel(md.stage));
            md.params.dI2R = replaceConstantByDistribution(md.params.dI2R, numel(md.stage));
            md.params.dI2H = replaceConstantByDistribution(md.params.dI2H, numel(md.stage));
            md.params.dH2R = replaceConstantByDistribution(md.params.dH2R, numel(md.stage));
            
            
            tGrids = md.params.tGrids(1) : md.dt : md.params.tGrids(end);
            md.record = [];
            for i = 1:numel(tGrids)
                info = updateStages(md);
                md.record = [md.record; info];
            end

            md.record.CumulativeI = cumsum(md.record.newI);
            md.record.CumulativeA = cumsum(md.record.newA);
            md.record.CumulativeR = cumsum(md.record.newR);
            md.record.CumulativeH = cumsum(md.record.newH);      
            md.record.tGrids = tGrids';

            function y = replaceConstantByDistribution(mean, N)
                x = rand(N,1);
                y = -log(1-x) * mean; % inverse transform of exponential distribution
                % y = round(y);
                % y(y==0) = 1;
            end
        end

    end %methods


end