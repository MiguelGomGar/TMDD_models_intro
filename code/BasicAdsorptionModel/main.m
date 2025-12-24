%% Comparison of Binding Isotherms: General vs. Langmuir
% Description: Visualizes the error introduced by the Langmuir 
% approximation when the drug concentration is not significantly higher 
% than receptor concentration.

clear all; clc; close all;

%% 1. Parameter Definition
R0 = 10;  % Total Receptor Concentration (e.g., nM)
KD = 2.0; % Equilibrium Dissociation Constant (e.g., nM)

% Create a range of Drug Concentrations (L0)
% We focus on the range close to R0 to see the TMDD effect clearly.
L0 = linspace(0, 50, 500); 

%% 2. Calculate Complex Concentrations (P)
% Model A: General Quadratic Solution (Exact) 
P_exact = binding_exact(L0, R0, KD);

% Model B: Standard Langmuir Approximation (Valid if L >>> R)
P_langmuir = binding_langmuir(L0, R0, KD);

%% 3. Visualization
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.5]);

% Plot Exact Model
plot(L0, P_exact, 'LineWidth', 2.5, 'Color', [0 0.447 0.741], ...
    'DisplayName', 'General Quadratic Model');
hold on;

% Plot Langmuir Approximation
plot(L0, P_langmuir, 'LineWidth', 2.5, 'DisplayName', 'Langmuir Approximation');

% Add reference line for R0 (Saturation level)
yline(R0, 'k:', 'Label', 'Saturation (R_{0})', ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

% Formatting
title('Ligand-Receptor Binding Isotherms', 'FontSize', 14);
xlabel('Total Drug Concentration ($L_0$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Complex Concentration ($P$)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 11);
grid on;
axis([0 max(L0) 0 R0*1.1]);