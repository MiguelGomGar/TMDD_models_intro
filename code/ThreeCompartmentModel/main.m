%% TMDD Three-Compartment Model Simulation (Depot Absorption)
% Description: Simulation of SC administration. Drug starts in the Depot 
% compartment and follows first-order absorption into the Central 
% compartment, where it undergoes TMDD and tissue distribution.

clear; clc; close all;

%% 1. Parameter Configuration
% Kinetic constants stored in a struct for clean passing
params.k_a   = 0.15;   % Absorption rate constant (Depot -> Central) [1/time]
params.k_eL  = 0.15;   % Elimination rate of free drug from Central [1/time]
params.k_eP  = 0.05;   % Elimination rate of the complex [1/time]
params.k_on  = 0.09;   % Association rate constant [1/(conc*time)]
params.k_off = 0.01;   % Dissociation rate constant [1/time]
params.k_in  = 1.0;    % Receptor synthesis rate [conc/time]
params.k_out = 0.2;    % Receptor degradation rate [1/time]
params.k_pt  = 0.08;   % Distribution: Plasma to Tissue [1/time]
params.k_tp  = 0.04;   % Distribution: Tissue to Plasma [1/time]

% Simulation Time
tspan = [0 80];

%% 2. Initial Conditions
% Dose is administered entirely into the Depot compartment.
Dose = 50.0; 

% Calculate baseline receptor concentration (Steady State)
R_ss = params.k_in / params.k_out;

% State Vector y0 = [Ld; Lc; Lt; R; P]
Ld0 = Dose;
Lc0 = 0.0;   
Lt0 = 0.0;
R0  = R_ss;
P0  = 0.0;

y0 = [Ld0; Lc0; Lt0; R0; P0];

%% 3. Numerical Integration using ode45
% Create anonymous function to pass parameters
ode_system = @(t, y) three_compartment_model(t, y, params);

% Solve the system
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y] = ode45(ode_system, tspan, y0, options);

% Extract results for visualization
Ld_res = y(:, 1); % Drug in Depot
Lc_res = y(:, 2); % Drug in Central
Lt_res = y(:, 3); % Drug in Tissue
R_res  = y(:, 4); % Free Receptor
P_res  = y(:, 5); % Complex

%% 4. Visualization of Results
figure('Name', 'Three-Compartment TMDD Model', 'Units', 'normalized', ...
    'Position', [0.1, 0.1, 0.8, 0.6]);

% --- Subplot 1: Pharmacokinetics (Drug Distribution) ---
subplot(1, 2, 1);
plot(t, Ld_res, 'LineWidth', 2, 'DisplayName', 'Depot (L_D)');
hold on;
plot(t, Lc_res, 'LineWidth', 2, 'DisplayName', 'Central Plasma (L_C)');
plot(t, Lt_res, 'LineWidth', 2, 'DisplayName', 'Peripheral Tissue (L_T)');
hold off;

title('Pharmacokinetics: Absorption & Distribution', 'FontSize', 12);
xlabel('Time', 'FontSize', 11);
ylabel('Free Drug Concentration', 'FontSize', 11);
legend('Location', 'best');
grid on;

% --- Subplot 2: Pharmacodynamics (Target Binding) ---
subplot(1, 2, 2);
plot(t, R_res, 'LineWidth', 2, 'DisplayName', 'Free Receptor (R)');
hold on;
plot(t, P_res, 'LineWidth', 2, 'DisplayName', 'Drug-Receptor Complex (P)');
yline(R0, ':', 'DisplayName', 'Baseline R_0');
hold off;

title('Pharmacodynamics: Target Engagement', 'FontSize', 12);
xlabel('Time', 'FontSize', 11);
ylabel('Concentration', 'FontSize', 11);
legend('Location', 'best');
grid on;