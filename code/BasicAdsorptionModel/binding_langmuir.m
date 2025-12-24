function P = binding_langmuir(L0, R0, KD)
    % BINDING_LANGMUIR Calculates complex concentration using the classic 
    % approximation.
    %
    % This function uses the standard Langmuir adsorption isotherm. It 
    % assumes that the concentration of the ligand greatly exceeds that of 
    % the receptor (L_free ~ L0).
    %
    % Arguments:
    %   L0 - Total drug concentration (vector or scalar).
    %   R0 - Total receptor concentration (scalar).
    %   KD - Equilibrium dissociation constant (scalar).
    %
    % Output:
    %   P  - Concentration of the drug-receptor complex.

    P = (R0 * L0) ./ (KD + L0);
end