%% Pharmacokinetic One-Compartment Model Simulation
% Description: Integration of the L-R-P system using ode45
clear; clc; close all;

%% 1. Parameter Configuration
% Kinetic constants
params.k_eL  = 0.15;  % Drug elimination rate
params.k_eP  = 0.05;  % Complex elimination rate
params.k_on  = 0.25;  % Binding rate
params.k_off = 0.08;  % Unbinding rate
params.k_in  = 1.20;  % Receptor production
params.k_out = 0.30;  % Receptor turnover

% Time span
tspan = [0 60];

%% 2. Initial Conditions
% Initial concentrations: [L0, R0, P0]
L0 = 12.0; % Initial bolus dose
R0 = 4.0;  % Initial receptor density (steady state)
P0 = 0.0;  % No complex formed at t = 0

y0 = [L0; R0; P0];

%% 3. Numerical Integration
% Define anonymous function to pass parameters
ode_system = @(t, y) one_compartment_model(t, y, ...
    params.k_eL, params.k_eP, params.k_on, params.k_off, params.k_in, params.k_out);

% Solver execution
[t, y] = ode45(ode_system, tspan, y0);

% Data arrangments
t = t';

L = y(:, 1)'; % free drug
R = y(:, 2)'; % free receptor
P = y(:, 3)'; % complex

%% 4. Results Visualization
figure('Units', 'normalized', 'Position', [0.1, 0.05, 0.7, 1]);

% Plotting the trajectories
plot(t, L, 'LineWidth', 2); 
hold on;
plot(t, R, 'LineWidth', 2);
plot(t, P, 'LineWidth', 2);

% Formatting
xlabel('Time', 'FontSize', 12);
ylabel('Concentration', 'FontSize', 12);
legend({'Free Drug (L)', 'Free Receptor (R)', 'Complex (P)'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);