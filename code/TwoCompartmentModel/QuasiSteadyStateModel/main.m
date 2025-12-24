%% TMDD Quasi-Steady-State (QSS) Model Simulation
clear; clc; close all;

%% 1. Parameters
p.k_eL  = 0.15;
p.k_eP  = 0.05;
p.k_in  = 1.0;
p.k_out = 0.2;
p.k_pt  = 0.08;
p.k_tp  = 0.04;
k_on    = 0.09;
k_off   = 0.01;

% Critical Difference: K_SS calculation
p.K_SS  = (k_off + p.k_eP) / k_on; 

%% 2. Initial Conditions & Run
y0 = [20.0; 0; p.k_in/p.k_out]; % [Ltot, Lt, Rtot]
tspan = [0 80];

[t, y] = ode45(@(t,y) model_QSS(t,y,p), tspan, y0);

%% 3. Post-processing (using K_SS)
L_tot = y(:,1); 
R_tot = y(:,3);
Lc = zeros(size(L_tot)); 
P = zeros(size(L_tot));

for i=1:length(t)
    term = L_tot(i) - R_tot(i) - p.K_SS;
    Lc(i) = 0.5 * (term + sqrt(term^2 + 4*p.K_SS*L_tot(i)));
    P(i) = (R_tot(i) * Lc(i)) / (p.K_SS + Lc(i));
end

%% 4. Visualization
figure('Name','QSS Model');
plot(t, Lc, 'LineWidth', 2, 'DisplayName', 'L_c (Free)'); 
hold on;
plot(t, P, 'LineWidth', 2, 'DisplayName', 'P (Complex)');
plot(t, R_tot, 'DisplayName', 'R_{tot}');
title('QSS Model Dynamics'); 
legend; 
grid on; 
xlabel('Time'); 
ylabel('Concentration');