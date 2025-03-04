function [err, coeff, misalign] = rls_function(lambda, filter_len, u, desired_sig, delta, fir_filter)
    %% Input arguments:
    % lambda = Forgetting factor, scalar (1x1)
    % filter_len = Filter length, scalar (number of coefficients in the adaptive filter)
    % u = Input signal, column vector (signal_len x 1)
    % desired_sig = Desired signal, column vector (signal_len x 1) (output of the unknown system + noise)
    % delta = Initial value for P(0) = delta^-1 * I, scalar (1x1)
    % fir_filter = True system impulse response, column vector (filter_len x 1)
    
    %% Output arguments:
    % err = A priori estimation error, column vector (signal_len x 1)
    % coeff = Final filter coefficients, column vector (filter_len x 1)
    % misalign = Misalignment in dB at each iteration, column vector (signal_len x 1)
    
    %% Initialization
    coeff = zeros(filter_len, 1); % Initialize adaptive filter coefficients to zero
    P = eye(filter_len) / delta;  % Initialize inverse correlation matrix with P(0) = delta^-1 * I
    
    % u and desired_sig are column vectors
    u = u(:);
    desired_sig = desired_sig(:);
    
    signal_len = length(u); % Length of the input signal
    
    err = zeros(signal_len, 1); % Initialize error vector to zero
    misalign = zeros(signal_len, 1); % Initialize misalignment vector to zero
    
    %% Main loop for RLS algorithm
    for n = filter_len:signal_len
        % Extract input vector: [u(n), u(n-1), ..., u(n-filter_len+1)]
        u_vec = u(n:-1:n-filter_len+1);
    
        % Compute Kalman gain vector (eq. 7 in the paper)
        kalman = (lambda^(-1) * P * u_vec) / (1 + lambda^(-1) * u_vec' * P * u_vec);
    
        % Compute a priori error (eq. 6 in the paper)
        err(n) = desired_sig(n) - coeff' * u_vec;
    
        % Update filter coefficients (eq. 8 in the paper)
        coeff = coeff + kalman * conj(err(n));
    
        % Update inverse correlation matrix (eq. 9 in the paper)
        P = (lambda^(-1) * P - lambda^(-1) * kalman * u_vec' * P);
    
        % Compute misalignment in dB
        misalign(n) = 20 * log10(norm(fir_filter - coeff)^2 / norm(fir_filter)^2);
    end
end