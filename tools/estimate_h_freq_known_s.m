function h_freq = estimate_h_freq_known_s(x,s,V,Cm)
threshold = 1e5;
max_error = 1e-4;

% Checking for numerical issues
if cond(diag(V'*s)) > 1e16
    disp('SKIP FOR NUMERICAL ISSUES!')
    h_freq = [];
    return
end

% Estimating h from analytical solution discard filters close to 0 
% to aboid numerical issues
h_freq = diag(V'*s)\V'*Cm'*x;
h_freq(abs(h_freq)>threshold) = threshold;

% Check if the optimization worked 
if norm(Cm'*x-V*diag(h_freq)*V'*s) > max_error
   h_freq = [];     
end


