clear; close all; clc;

%% Classic Model
syms beta c q kappa p omega omegap alpha gamma gammap gammapp tau;
syms S E I A R H N;
var = [E, I, A, H];

f = [beta * S * (I + kappa*A) / N; 0; 0; 0];
v = [p*omega*E + (1-p)*omegap*E;...
    -(1-p)*omegap*E + alpha*gamma*I + (1-alpha)*gammap*I;...
    -p*omega*E + gammapp*A;...
    -(1-alpha)*gammap*I + tau*H];

F = jacobian(f, var);
V = jacobian(v, var);

M = F * inv(V);

[~, D] = eig(M);
D = diag(D);
%Reff = simplify(D(end));
Reff = D(end);
tex = latex(Reff)

%% Improved Model
syms beta c q kappa p omega omegap alpha gamma gammap gammapp tau;
syms S E1 E2 I1 I2 A R H;
var2 = [E1, E2, I1, I2, A, H];

f2 = [p*beta*S*(I1 + I2 + kappa*A) / N;...
    (1-p)*beta*S*(I1 + I2 + kappa*A) / N;...
    0; 0; 0; 0];
v2 = [omega*E1;...
    omegap*E2;...
    -alpha*omegap*E2 + gamma*I1;...
    -(1-alpha)*omegap*E2 + gammap*I2;...
    -omega*E1 + gammapp*A;...
    -gammap*I2 - tau*H];

F2 = jacobian(f2, var2);
V2 = jacobian(v2, var2);

M2 = F2 * inv(V2);

[~, D2] = eig(M2);
D2 = diag(D2);
Reff2 = simplify(D2(end));
%Reff2 = D2(end);
tex2 = latex(Reff2)


%% Difference
simplify(Reff - Reff2)