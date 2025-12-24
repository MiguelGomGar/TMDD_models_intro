function dy_dt = three_compartment_model(t, y, params)
    % THREE_COMPARTMENT_MODEL TMDD with Tissue Distribution and Depot 
    % Absorption.
    %
    % This function extends the two-compartment model by adding a depot 
    % compartment for extravascular administration (e.g., Subcutaneous).
    %
    % The drug flows: Depot -> Central <-> Tissue.
    % Binding occurs ONLY in the Central compartment.
    %
    % Arguments:
    %   t      - Time vector.
    %   y      - State vector [Ld; Lc; Lt; R; P]
    %            y(1) Ld: Drug in Depot compartment (New).
    %            y(2) Lc: Free drug in Central compartment.
    %            y(3) Lt: Free drug in Tissue compartment.
    %            y(4) R:  Free receptor concentration.
    %            y(5) P:  Drug-receptor complex concentration.
    %   params - Struct containing kinetic constants, including 'k_a'.

    % 1. Map state variables
    Ld = y(1); % Depot
    Lc = y(2); % Central
    Lt = y(3); % Tissue
    R  = y(4); % Receptor
    P  = y(5); % Complex

    % 2. Extract parameters
    k_a   = params.k_a;   % Absorption rate constant (New)
    k_eL  = params.k_eL;
    k_eP  = params.k_eP;
    k_on  = params.k_on;
    k_off = params.k_off;
    k_in  = params.k_in;
    k_out = params.k_out;
    k_pt  = params.k_pt;  % Plasma to Tissue
    k_tp  = params.k_tp;  % Tissue to Plasma

    % 3. Differential Equations
    
    % Depot Absorption
    dLd_dt = -k_a * Ld;

    % dLc/dt: Central (Input from Depot - Elim - Binding - Dist)
    dLc_dt = k_a*Ld - (k_eL + k_pt)*Lc - k_on*Lc*R + k_off*P + k_tp*Lt;

    % dLt/dt: Tissue Distribution
    dLt_dt = k_pt*Lc - k_tp*Lt;

    % dR/dt: Receptor Turnover & Binding
    dR_dt = k_in - k_out*R - k_on*Lc*R + k_off*P;

    % dP/dt: Complex Dynamics
    dP_dt = k_on*Lc*R - k_off*P - k_eP*P;

    % 4. Return derivatives
    dy_dt = [dLd_dt; dLc_dt; dLt_dt; dR_dt; dP_dt];
end