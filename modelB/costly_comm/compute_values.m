function [lambda, D_k1, N_k1] = compute_values(k, beta, var,a)

    % We store all the tolerances in this function rather than passing them
    % around
    AbsTol = 1e-8;
    RelTol = 1e-6;

    % The epsilon is used to calculate the derivative numerically
    %for cauchy distribution
    scaling_par=1;
    location_par=0;
    
    epsilon = 1e-4;

    RHS_M  = @(e) e.^0;
    RHS_L  = @(e) e.^2;
    %kernel = @(e,x) (1/sqrt(2*pi*var)) * exp(-(x-a*e).^2 / (2*var)); %Gaussian
    %kernel=@(e,x) (1/pi)*(scaling_par./((x-e-location_par).^2 + (scaling_par)^2));  %Cauchy

    %For uniform distribution
    state_max=15;

    kernel=@(e,x) (abs((x-e)) < state_max)*(1/(2*state_max)); % Uniform


    k1 = k;
    k2 = k1 + epsilon;

    if k1==0
        L_k1 = 0;
        M_k1 = 1;
        D_k1 = 0;
        N_k1 = 1;
    else
        M_k1 = Fie0(beta, k1, kernel, RHS_M, AbsTol, RelTol);
        L_k1 = Fie0(beta, k1, kernel, RHS_L, AbsTol, RelTol);
        
        D_k1 = L_k1 / M_k1;
        N_k1 =    1 / M_k1 - (1-beta);
    end

    M_k2 = Fie0(beta, k2, kernel, RHS_M, AbsTol, RelTol);
    L_k2 = Fie0(beta, k2, kernel, RHS_L, AbsTol, RelTol);
    
    D_k2 = L_k2 / M_k2;
    N_k2 =    1 / M_k2 - (1-beta);

    lambda = -(D_k2 - D_k1) / (N_k2 - N_k1);

    % An alternative method to compute lambda is as follows:
    %    
    % L_dot = (L_k2 - L_k1) / epsilon;
    % M_dot = (M_k2 - M_k1) / epsilon;
    %
    % D_dot = (M_k1*L_dot - L_k1*M_dot) / ((M_k1)^2);
    % N_dot = -(M_dot/(M_k1)^2);
    % 
    % lambda = -(D_dot/N_dot);
    
end 
