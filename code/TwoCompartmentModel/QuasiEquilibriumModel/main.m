%% TMDD Quasi-Equilibrium (QE) Model Simulation
clear; clc; close all;

%% 1. Parameters
p.k_eL  = 0.15;
p.k_eP  = 0.05;
p.k_in  = 1.0;
p.k_out = 0.2;
p.k_pt  = 0.08;
p.k_tp  = 0.04;

% Derived parameters for QE
k_on  = 0.09;
k_off = 0.01;
p.KD  = k_off / k_on; % Equilibrium constant

%% 2. Initial Conditions
L0 = 20.0; % Initial Bolus
R0 = p.k_in / p.k_out;

% Initial State: [L_tot, L_t, R_tot]
% At t = 0, P=0 implies L_tot = Lc0 and R_tot = R0
y0 = [L0; 0; R0];

tspan = [0 80];

%% 3. Simulation
[t, y] = ode45(@(t,y) model_QE(t,y,p), tspan, y0);

L_tot_res = y(:,1);
L_t_res   = y(:,2);
R_tot_res = y(:,3);

%% 4. Post-processing: Recover Free Concentrations (Lc, R, P)
Lc_res = zeros(size(L_tot_res));
P_res  = zeros(size(L_tot_res));
R_res  = zeros(size(L_tot_res));

for i = 1:length(t)
    L_tot = L_tot_res(i);
    R_tot = R_tot_res(i);
    term = L_tot - R_tot - p.KD;
    Lc_res(i) = 0.5 * (term + sqrt(term^2 + 4*p.KD*L_tot));
    P_res(i)  = L_tot - Lc_res(i);
    R_res(i)  = R_tot - P_res(i);
end

%% 5. Plotting
figure('Name','QE Model');
subplot(1,2,1);
plot(t, L_tot_res, 'DisplayName', 'L_{tot}'); 
hold on;
plot(t, Lc_res, 'DisplayName', 'L_c (Free)');
title('Drug Kinetics (QE)'); 
legend; 
grid on;

subplot(1,2,2);
plot(t, R_tot_res, 'DisplayName', 'R_{tot}'); 
hold on;
plot(t, R_res, 'DisplayName', 'R (Free)');
plot(t, P_res, 'DisplayName', 'P (Complex)');
title('Target Dynamics (QE)'); 
legend; 
grid on;