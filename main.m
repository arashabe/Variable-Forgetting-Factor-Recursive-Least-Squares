clear;
close all;
clc;

%% Parameters
filter_len = 64; % Filter length (as per the paper)
iterations = 15000; % Number of iterations
u = wgn(iterations, 1, 0); % Input signal: White Gaussian Noise (generates an m-by-n matrix of white Gaussian noise samples in volts)
signal_len = length(u); % Input signal length
frequency = 0.5; % Frequency of fundamental tone (0 < freq < 1)
fir_filter = fir1(filter_len-1, frequency)'; % Uses a Hamming window to design a 64th-order lowpass FIR filter (https://it.mathworks.com/help/signal/ref/fir1.html)
SNR = 20; % Signal-to-noise ratio of the tone
delta = 1; % Initial value for P(0) = delta^-1 * I, scalar (1x1)

%% Main loop for two types of input signals
for i = 1:2
    % For the second iteration, change the input signal to an AR(1) process
    if i == 2
        u = filter([1], [1 -0.9], u); % Input signal is AR(1) Process
    end
    
    % Compute the output of the unknown system
    y = zeros(signal_len, 1); % Initialize output vector
    for n = filter_len:signal_len
        u_vec = u(n:-1:n-filter_len+1); % Input signal vector X(n)
        y(n) = fir_filter' * u_vec; % Output of the unknown system
    end
    
    % Add noise to the desired signal to achieve the specified SNR
    desired_sig = awgn(y, SNR); % Desired signal: output of the unknown system corrupted by white Gaussian noise with 20 dB SNR
    
    %% RLS Algorithm
    lambda = 1 - 1 / (3 * filter_len); % Forgetting factor for RLS
    [err_rls, coeff_rls, misalign_rls] = rls_function(lambda, filter_len, u, desired_sig, delta, fir_filter); % Run RLS algorithm
    
    %% VFF-RLS Algorithm
    Ka = 2; % Parameter for exponential window
    Kb = 5 * Ka; % Parameter for exponential window
    [err_vff, coeff_vff, misalign_vff, lambda_evol, cond_num] = vff_rls_function(filter_len, u, desired_sig, delta, fir_filter, Ka, Kb); % Run VFF-RLS algorithm
    
    %% Plot results for White Gaussian Noise input
    if i == 1
        figure;
        plot(filter_len:signal_len, misalign_rls(filter_len:signal_len), filter_len:signal_len, misalign_vff(filter_len:signal_len));
        xlabel('Iterations');
        ylabel('Misalignment (dB)');
        legend('RLS', 'VFF-RLS');
        title('Input Signal is White Gaussian Noise');
        
        figure;
        plot(filter_len:signal_len, lambda_evol(filter_len:signal_len));
        xlabel('Iterations');
        ylabel('\lambda');
        legend('VFF-RLS');
        title('Evolution of Forgetting Factor - White Gaussian Noise');
    end
    
    %% Plot results for AR(1) Process input
    if i == 2
        figure;
        plot(filter_len:signal_len, misalign_rls(filter_len:signal_len), filter_len:signal_len, misalign_vff(filter_len:signal_len));
        xlabel('Iterations');
        ylabel('Misalignment (dB)');
        legend('RLS', 'VFF-RLS');
        title('Input Signal is AR(1) Process');
        
        figure;
        plot(filter_len:signal_len, lambda_evol(filter_len:signal_len));
        xlabel('Iterations');
        ylabel('\lambda');
        legend('VFF-RLS');
        title('Evolution of Forgetting Factor - AR(1) Process');
    end
end