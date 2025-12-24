function dy_dt = model_QSS(t, y, p)
    % MODEL_QSS Quasi-Steady-State TMDD approximation.
    %
    % Similar to QE but uses K_SS instead of K_D.
    % Valid when k_eP is significant compared to k_off.
    %    
    % The equations structure is identical to QE, just swapping constants.
    % We reuse the logic but p.K_SS must be provided.
    
    L_tot = y(1);
    L_t   = y(2);
    R_tot = y(3);

    Kss = p.K_SS; % The only difference

    % Quadratic solution for Free Drug (Lc) using Kss
    term_diff = L_tot - R_tot - Kss;
    discriminant = term_diff^2 + 4 * Kss * L_tot;
    Lc = 0.5 * (term_diff + sqrt(discriminant));

    % Derivatives
    dLt_dt = p.k_pt * Lc - p.k_tp * L_t;

    complex_term = (R_tot * Lc) / (Kss + Lc); % Proportional to P

    dLtot_dt = -(p.k_eL + p.k_pt)*Lc - p.k_eP * complex_term + p.k_tp*L_t;
    dRtot_dt = p.k_in - p.k_out*R_tot - (p.k_eP - p.k_out) * complex_term;

    dy_dt = [dLtot_dt; dLt_dt; dRtot_dt];
end