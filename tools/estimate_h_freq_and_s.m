function [h_freq,s,error] = estimate_h_freq_and_s(x,V,S)
% Estimate the frequency response of the filter coefficients h and the
% sparse signal s by solving the problem min(norm(H_inv*x-s)^2+norm(s,1))
% using alternating minimization and extracting the freq. response of h
% from H_inv, assuming H = inv(H_inv) and H=V*diag(h_freq)*V'
% S represent the number of seeds seeked during the minimization.

n_iters = 200;
lambda = 0.1;
threshold = 1e-5;
max_error = 1e-4;

N = length(x);
s = zeros(N,1);
h_inv = ones(N,1)*0.9;
error = zeros(n_iters,1);

% Checking for numerical issues
if cond(diag(V'*x)) > 1e15
    h_freq = [];
    return
end

for k=1:n_iters
    % Estimating s* from analytical solution
    for i=1:N
        e_i = [zeros(i-1,1); 1; zeros(N-i,1)];
        alpha_i = e_i'*V*diag(V'*x)*h_inv;
        if abs(alpha_i) > lambda
            s(i) = alpha_i - lambda*sign(s(i));
        else
            s(i) = 0;
        end
    end

    % Estimating h_inv from analytical solution discard filters close to 0 
    % to aboid numerical issues
    h_inv = diag(V'*x)\V'*s;
    h_inv(abs(h_inv)<threshold) = threshold;
    error(k) = norm(V*diag(V'*x)*h_inv-s)^2+norm(s,1);

    % If s* is sparse enough stop algorithm
    if sum(abs(s)>0) <= S
        break
    end
end

% Check if the optimization worked 
h_freq = h_inv.^-1;
if norm(x-V*diag(h_freq)*V'*s) > max_error
    h_freq = [];     
end