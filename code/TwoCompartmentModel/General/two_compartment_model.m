function dy_dt = two_compartment_model(t, y, params)
    % TWO_COMPARTMENT_MODEL Calculates derivatives for TMDD with tissue 
    % distribution.
    %
    % This function models the kinetics where the drug distributes between 
    % a central compartment and a peripheral tissue compartment. Binding to
    % the target receptor occurs ONLY in the central compartment.
    %
    % Arguments:
    %   t      - Time vector.
    %   y      - State vector [Lc; Lt; R; P]
    %            y(1) Lc: Free drug concentration in Central compartment.
    %            y(2) Lt: Free drug concentration in Tissue compartment.
    %            y(3) R:  Free receptor concentration.
    %            y(4) P:  Drug-receptor complex concentration.
    %   params - Struct containing kinetic constants:
    %            .k_eL  (Elimination of free drug from central)
    %            .k_eP  (Elimination of complex)
    %            .k_on  (Association rate)
    %            .k_off (Dissociation rate)
    %            .k_in  (Receptor synthesis)
    %            .k_out (Receptor degradation)
    %            .k_pt  (Transfer: Plasma -> Tissue)
    %            .k_tp  (Transfer: Tissue -> Plasma)
    
    % 1. Map state variables for readability
    Lc = y(1); % Drug in Central Cpt
    Lt = y(2); % Drug in Tissue Cpt
    R  = y(3); % Free Receptor
    P  = y(4); % Complex

    % 2. Extract parameters
    k_eL  = params.k_eL;
    k_eP  = params.k_eP;
    k_on  = params.k_on;
    k_off = params.k_off;
    k_in  = params.k_in;
    k_out = params.k_out;
    k_pt  = params.k_pt;
    k_tp  = params.k_tp;

    % 3. Differential Equations
    
    % dLc/dt: Central Drug (Elimination, Binding, Distribution)
    dLc_dt = -k_eL*Lc - k_on*Lc*R + k_off*P - k_pt*Lc + k_tp*Lt;

    % dLt/dt: Tissue Drug (Distribution in/out)
    dLt_dt = k_pt*Lc - k_tp*Lt;

    % dR/dt: Free Receptor (Turnover + Binding)
    dR_dt = k_in - k_out*R - k_on*Lc*R + k_off*P;

    % dP/dt: Complex (Formation - Dissociation - Elimination)
    dP_dt = k_on*Lc*R - k_off*P - k_eP*P;

    % 4. Return derivatives column vector
    dy_dt = [dLc_dt; dLt_dt; dR_dt; dP_dt];
end