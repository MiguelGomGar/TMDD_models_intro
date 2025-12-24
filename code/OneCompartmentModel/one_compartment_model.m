function dy_dt = one_compartment_model(t, y, k_eL, k_eP, k_on, k_off, k_in, k_out)
    
    % Calculates the derivative of the concentration of free drug,
    % available receptors and drug-receptor complex with respect to time
    % given a a single bolus infusion of the drug into the central 
    % compartment.

    % Parameters:
    %
    %   t: time vector.
    %   
    %   y: vector that wrappes up the variables L, R, P.
    %
    %       L: concentration of the free drug.
    %
    %       R: concentration of free receptor.
    %
    %       P: concentration of the drug-receptor complex.
    %
    %   k_eL: kinetic constant for the free drug elimination process.
    % 
    %   k_eP: kinetic constant for the product elimination process.
    %
    %   k_on: kinetic constant for the drug-complex association process.
    %
    %   k_off: kinetic constant for the drug-complex dissociation process.
    %
    %   k_in: kinetic constant for the receptor creation process.
    %
    %   k_out: kinetic constant for the receptor elimination process.
    %
    %
    % Outputs:
    %
    %   dL_dt: derivative of the concentration of free drug with respect to
    %   time
    %
    %   dR_dt: derivative of the concentration of free receptor with 
    %   respect to time.
    %
    %   dP_dt: derivative of the concentration of the drug-receptor complex
    %   with respect to time.
    
    L = y(1);
    R = y(2);
    P = y(3);

    dL_dt = -k_eL*L - k_on*L*R + k_off*P;
    dR_dt = k_in - k_out*R - k_on*L*R + k_off*P;
    dP_dt = k_on*L*R - k_off*P -k_eP*P;

    dy_dt = [dL_dt; dR_dt; dP_dt];

end