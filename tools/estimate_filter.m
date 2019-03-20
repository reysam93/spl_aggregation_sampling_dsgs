% Using the loaded data tries to estimate the filter coefficients and the
% sparse signal s so x can be represented as a diffusion process of the
% form x=Hs. The data may be saved so it can be used fot the recontruction
% algorithms.
load('..\dataset_additional\graphs.mat');

save_data = false;
gso = 'L'; % Options A or L
n_graphs = length(cell_X);
max_nodes = 75;
min_nodes = 20;
attr = 6;
S=3;

tic
hh_freq = cell(n_graphs,1);
estimated_ss = cell(n_graphs,1);
vv = cell(n_graphs,1);
ee = cell(n_graphs,1);
err_no_skip = zeros(n_graphs,1);
mean_err = zeros(n_graphs,1);
mean_sparsity = zeros(n_graphs,1);
skipped = 0;
for i=1:n_graphs
    % Select data from graph i
    x = cell_X{i}(:,attr);
    A = cell_A{i};
    N = length(A);
    hh_freq{i} = [];
    
    % Eigendecomposition of the GSO
    if strcmp(gso,'L')
        G = gsp_graph(A);
        G = gsp_compute_fourier_basis(G);
        e = G.e;
        V = G.U;
    elseif strcmp(gso,'A')
        [V,D] = eig(A);
        e = diag(D);
    end
    
    % Limiting graph size
    if N < min_nodes || N > max_nodes
        skipped = skipped + 1;
        continue
    end
    
    % Estimating freq. response of h and sparse s
    [h_freq,s_est]=estimate_h_freq_and_s(x,V,S);

    if isempty(h_freq)
        skipped = skipped + 1;
        continue
    end
    estimated_ss{i} = s_est;
    hh_freq{i} = h_freq;
    vv{i} = V;
    ee{i} = e;
    err_no_skip(i) = norm(V*diag(h_freq)*V'*s_est-x);
    disp(['Graph ' num2str(i) '   N: ' num2str(length(A)) '   Sparsity: '... 
        num2str(sum(abs(s_est)>0)) '   MSE: ' num2str(err_no_skip(i))])
end

% Remove skipped graphs from data
SS = cell(n_graphs-skipped,1);
XX = cell(n_graphs-skipped,1);
HH_Freq = cell(n_graphs-skipped,1);
VV = cell(n_graphs-skipped,1);
EE = cell(n_graphs-skipped,1);
err_skip = zeros(n_graphs-skipped,1);
max_h = zeros(n_graphs-skipped,1);
sparsity = zeros(n_graphs-skipped,1);
cont = 1;
for j=1:n_graphs
    if isempty(hh_freq{j})
        continue
    end
    XX{cont} = cell_X{j}(:,attr);
    SS{cont} = estimated_ss{j};
    HH_Freq{cont} = hh_freq{j};
    VV{cont} = vv{j};
    EE{cont} = ee{j};
    max_h(cont) = max(hh_freq{j});
    err_skip(cont) = err_no_skip(j);
    sparsity(cont) = sum(abs(estimated_ss{j})>0);
    cont = cont + 1;
end
toc

disp(['Mean sparsity:   ' num2str(mean(sparsity))])
disp(['Mean MSE:        ' num2str(mean(err_skip))])
disp(['Selected graphs: ' num2str(length(XX))])
disp(['Mean Max h       ' num2str(mean(max_h))])

if save_data
    save(['..\dataset\all_data_eig' gso '_x' num2str(attr) '_S' num2str(S)],...
            'EE', 'VV', 'XX', 'SS', 'HH_Freq')
end