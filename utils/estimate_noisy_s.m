function s = estimate_noisy_s(y_m, Cm, Theta, Rw_m, gamma, verbose)
% Estimates the original signal s from the observed (sampled) signal y_mi
% with noise.
% Arguments:
%   - y_m: sampled signal
%   - Cm: sampling matrix
%   - Theta: observation matrix
%   - Rw_m: sampled covariance matrix of the noise
%   - gamma: regularizing parameter for the norm one

if ~exist('verbose', 'var')
    verbose = false;
end

N = size(Cm,2);

if verbose
    disp('Begining optimization')
end

cvx_begin quiet
    variable s(N)
    minimize (norm(Rw_m^(-1/2)*(y_m-Cm*Theta*s)) + gamma*norm(s,1))
cvx_end

if verbose
    fprintf('opt val: %d\n', cvx_optval)
    fprintf('\tnorm2: %d\n', norm(Rw_m^(-1/2)*(y_m-Cm*Theta*s)))
    fprintf('\tnorm1: %d\n', norm(s,1)*gamma)
end
