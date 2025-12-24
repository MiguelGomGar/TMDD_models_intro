function dy_dt = model_QE(t, y, p)
    % MODEL_QE Quasi-Equilibrium TMDD approximation.
    %
    % Assumes instantaneous equilibrium between L, R, and P.
    % Solves for TOTAL concentrations (L_tot, R_tot) instead of free ones.
    %
    % Arguments:
    %   t - Time vector
    %   y - State vector: [L_tot; L_t; R_tot]
    %       L_tot: Total drug in central cpt (Lc + P)
    %       L_t:   Drug in tissue cpt
    %       R_tot: Total receptor concentration (R + P)
    %   p - Parameters struct (needs KD, k_eL, k_eP, k_pt, k_tp, k_in, k_out)

    % 1. Unpack states
    L_tot = y(1);
    L_t   = y(2);
    R_tot = y(3);

    % 2. Algebraic step: Calculate free Lc using quadratic equation
    KD = p.KD;
    term_diff = L_tot - R_tot - KD;
    discriminant = term_diff^2 + 4 * KD * L_tot;
    Lc = 0.5 * (term_diff + sqrt(discriminant));

    % 3. Differential Equations
    
    % dL_t/dt: Tissue distribution
    dLt_dt = p.k_pt * Lc - p.k_tp * L_t;

    % Common term: P * k_eP approximation
    % Since P = (R_tot * Lc) / (KD + Lc), the elimination term is:
    complex_elim = (R_tot * p.k_eP * Lc) / (KD + Lc);

    % dL_tot/dt: Total Central Drug
    dLtot_dt = -(p.k_eL + p.k_pt)*Lc - complex_elim + p.k_tp*L_t;

    % dR_tot/dt: Total Receptor
    dRtot_dt = p.k_in - p.k_out*R_tot - (p.k_eP - p.k_out) * (R_tot * Lc)/(KD + Lc);

    dy_dt = [dLtot_dt; dLt_dt; dRtot_dt];
end