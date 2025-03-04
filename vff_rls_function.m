function [err, coeff, misalign, lambda_evol, cond_num] = vff_rls_function(filter_len, u, desired_sig, delta, fir_filter, Ka, Kb)
    %% Input arguments:
    % filter_len = Filter length, scalar (number of coefficients in the adaptive filter)
    % u = Input signal, column vector (signal_len x 1)
    % desired_sig = Desired signal, column vector (signal_len x 1) (output of the unknown system + noise)
    % delta = Initial value for P(0) = delta^-1 * I, scalar (1x1)
    % fir_filter = True system impulse response, column vector (filter_len x 1)
    % Ka, Kb = Parameters for exponential windows, scalars (1x1)
    
    %% Output arguments:
    % err = A priori estimation error, column vector (signal_len x 1)
    % coeff = Final filter coefficients, column vector (filter_len x 1)
    % misalign = Misalignment in dB at each iteration, column vector (signal_len x 1)
    % lambda_evol = Evolution of the forgetting factor, column vector (signal_len x 1)
    % cond_num = Condition number of the input signal covariance matrix, column vector (signal_len x 1)
    
    %% Compute exponential window parameters
    alpha = 1 - 1 / (Ka * filter_len); % Weighting factor for power estimates (eq. 14 and 15)
    beta = 1 - 1 / (Kb * filter_len);  % Weighting factor for noise power estimate (eq. 16)
    
    %% u and desired_sig are column vectors
    u = u(:);
    desired_sig = desired_sig(:);
    
    %% Initialize variables
    signal_len = length(u); % Length of the input signal
    coeff = zeros(filter_len, 1); % Initialize adaptive filter coefficients to zero
    P = eye(filter_len) / delta;  % Initialize inverse correlation matrix with P(0) = delta^-1 * I
    err = zeros(signal_len, 1);   % Initialize error vector to zero
    sigma_e = 1;         % Initialize power estimate of a priori error (sigma_e^2)
    sigma_q = 1;         % Initialize power estimate of q (sigma_q^2)
    sigma_v = 1;         % Initialize power estimate of system noise (sigma_v^2)
    lambda = 1;          % Initialize forgetting factor
    R = zeros(filter_len, filter_len); % Initialize input signal covariance matrix
    lambda_evol = zeros(signal_len, 1); % Initialize vector to store forgetting factor evolution
    cond_num = zeros(signal_len, 1);    % Initialize vector to store condition number evolution
    
    %% Main loop for VFF-RLS algorithm
    for n = filter_len:signal_len
        % Extract input vector: [u(n), u(n-1), ..., u(n-filter_len+1)]
        u_vec = u(n:-1:n-filter_len+1);
    
        % Compute Kalman gain vector (eq. 7)
        kalman = (lambda^(-1) * P * u_vec) / (1 + lambda^(-1) * u_vec' * P * u_vec);
    
        % Compute a priori error (eq. 6)
        err(n) = desired_sig(n) - coeff' * u_vec;
    
        % Update filter coefficients (eq. 8)
        coeff = coeff + kalman * conj(err(n));
    
        % Intermediate variable q (used in eq. 15)
        q = u_vec' * P * u_vec;
    
        % Update inverse correlation matrix (eq. 9)
        P = lambda^(-1) * P - lambda^(-1) * kalman * u_vec' * P;
    
        % Update power estimates
        sigma_e = alpha * sigma_e + (1 - alpha) * (err(n)^2); % Update sigma_e^2 (eq. 14)
        sigma_q = alpha * sigma_q + (1 - alpha) * (q^2);      % Update sigma_q^2 (eq. 15)
        sigma_v = beta * sigma_v + (1 - beta) * (err(n)^2);   % Update sigma_v^2 (eq. 16)
    
        % Compute gamma ratio (used to control the forgetting factor)
        gamma = sqrt(sigma_e) / sqrt(sigma_v);
    
        % Update forgetting factor based on gamma (eq. 17 and 18)
        if (gamma > 1) && (gamma <= 2)
            lambda = 1; % Lambda is 1 if gamma is in the specified range (eq. 17)
        else
            % Otherwise, update lambda using the formula from eq. 18
            lambda = min((sqrt(sigma_q) * sqrt(sigma_v)) / (1e-8 + abs(sqrt(sigma_e) - sqrt(sigma_v))), 1);
        end
    
        % Compute misalignment in dB
        misalign(n) = 20 * log10(norm(fir_filter - coeff)^2 / norm(fir_filter)^2);
    
        % Update input signal covariance matrix
        R = R + (lambda^(n-filter_len)) * (u_vec * u_vec');
    
        % Compute condition number of the input signal covariance matrix
        cond_num(n) = cond(R, 2);
    
        % Store the current value of lambda
        lambda_evol(n) = lambda;
    end
end