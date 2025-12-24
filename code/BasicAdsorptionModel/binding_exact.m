function P = binding_exact(L0, R0, KD)
    % BINDING_EXACT Calculates the equilibrium concentration of the complex
    % (P).
    %
    % This function uses the quadratic solution which represents the exact 
    % binding equilibrium without assuming that the concentration of the 
    % free drug is much larger than the concentration of receptors. It is 
    % valid for high-affinity compounds exhibiting TMDD.
    %
    % Arguments:
    %   L0 - Total initial free drug concentration (vector or scalar).
    %   R0 - Total receptor concentration (scalar).
    %   KD - Equilibrium dissociation constant (scalar).
    %
    % Output:
    %   P  - Concentration of the drug-receptor complex at equilibrium.

    term_sum = KD + L0 + R0;
    discriminant = term_sum.^2 - 4 * L0 * R0;
    
    P = 0.5 * (term_sum - sqrt(discriminant));
end