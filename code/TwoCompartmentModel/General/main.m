%% TMDD Two-Compartment Model Simulation
% Description: Integration of the full 2-compartment system (Binding in 
% Central).

clear; clc; close all;

%% 1. Parameter Configuration
% Defined in struct 'p' for cleaner passing to functions
p.k_eL  = 0.15;   % Central elimination of drug
p.k_eP  = 0.05;   % Complex elimination
p.k_on  = 0.09;   % Association rate (slower to see distribution effects)
p.k_off = 0.01;   % Dissociation rate
p.k_in  = 1.0;    % Receptor synthesis rate
p.k_out = 0.2;    % Receptor degradation rate

% Distribution constants (Plasma <-> Tissue)
p.k_pt  = 0.08;   % Transfer: Plasma to Tissue
p.k_tp  = 0.04;   % Transfer: Tissue to Plasma

% Time span
tspan = [0 100];

%% 2. Initial Conditions
Lc0 = 20.0;             % Initial Bolus in Central Compartment
Lt0 = 0.0;              % Initially 0 in Tissue
R0  = p.k_in / p.k_out; % Steady state receptor density
P0  = 0.0;              % No complex at t = 0

% State vector: [Lc, Lt, R, P]
y0 = [Lc0; Lt0; R0; P0];

%% 3. Numerical Integration
% Anonymous function wrapper to pass the 'p' struct
ode_system = @(t, y) two_compartment_model(t, y, p);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y] = ode45(ode_system, tspan, y0, options);

% Extract results for clarity
Lc_res = y(:, 1);
Lt_res = y(:, 2);
R_res  = y(:, 3);
P_res  = y(:, 4);

%% 4. Visualization

% Plot 1: Drug Distribution (Central vs Tissue)
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.5]);
plot(t, Lc_res, 'LineWidth', 2, 'DisplayName', 'Central Drug (L_C)');
hold on;
plot(t, Lt_res, 'LineWidth', 2, 'DisplayName', 'Tissue Drug (L_T)');
xlim([0, 50]);
xlabel('Time');
ylabel('Concentration');
title('Drug Distribution Kinetics');
legend('Location', 'bestoutside');
grid on;
set(gca, 'FontSize', 11);

% Plot 2: Receptor & Complex Dynamics
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.5]);
plot(t, R_res, 'LineWidth', 2, 'DisplayName', 'Free Receptor (R)');
hold on;
plot(t, P_res, 'LineWidth', 2, 'DisplayName', 'Complex (P)');
yline(R0, ':', 'Baseline R_0', 'HandleVisibility', 'off'); % Reference line
xlabel('Time');
ylabel('Concentration');
title('Target Binding Dynamics');
legend('Location', 'bestoutside');
grid on;
set(gca, 'FontSize', 11);