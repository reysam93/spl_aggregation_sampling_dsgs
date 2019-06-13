addpath('./utils')

% Load data from dataset
load('dataset\all_data_eigL_x6');

% Control variables
percentile = 25;
threshold = 5e-4;
gamma = 1e-2;
MM = 7:15;
n_graphs = min(100,length(XX));

g_prctile_mse = zeros(n_graphs,length(MM),4);
tic
for i=1:n_graphs
    % Select data from graph i
    h_freq = HH_Freq{i};
    e = EE{i};
    V = VV{i};
    x = XX{i};
    s = SS{i};
    N = length(h_freq);
    Psi = fliplr(vander(e));
    H = V*diag(h_freq)*V';   
    disp(['Graph: ' num2str(i) '   N_nodes: ' num2str(N)])
    
    mse_x_dsgs_1 = zeros(N,length(MM));
    mse_x_dsgs_2 = zeros(N,length(MM));
    mse_x_dsgs_5 = zeros(N,length(MM));
    mse_bgs = zeros(N,length(MM));
    for j=1:length(MM)
        M = MM(j);
        fprintf('\tM: %d\n', M)
        Cm = eye(M,N);
        Cms = rand_sampling_matrices(N,M,N);
        parfor node1=1:N
            % Random selection of samping nodes
            nodes_aux = 1:N;
            nodes_aux(node1) = [];
            nodes = [node1 nodes_aux(randperm(length(nodes_aux)-1,M-1))];
            
            % Assuming DSGS - 1 sampling node
            Xi_i = Psi'*diag(V(node1,:)'.*h_freq)*V';
            y_mi = Cm*Xi_i*s;
            s_dsgs_1 = estimate_noisy_s(y_mi,Cm,Xi_i,eye(M),gamma);
            mse_x_dsgs_1(node1,j) = norm(x-H*s_dsgs_1)^2;
            
            % Assuming DSGS - 2 sampling nodes
            [Cm_eq_2,Xi_eq_2] = eq_samp_matrices_k_nodes(N,M,nodes(1:2),Psi,V,h_freq);
            y_mi = Cm_eq_2*Xi_eq_2*s;
            s_dsgs_2 = estimate_noisy_s(y_mi,Cm_eq_2,Xi_eq_2,eye(M),gamma);
            mse_x_dsgs_2(node1,j) = norm(x-H*s_dsgs_2)^2;
            
            % Assuming DSGS - 5 sampling nodes
            [Cm_eq_5,Xi_eq_5] = eq_samp_matrices_k_nodes(N,M,nodes(1:5),Psi,V,h_freq);
            y_mi = Cm_eq_5*Xi_eq_5*s;
            s_dsgs_5 = estimate_noisy_s(y_mi,Cm_eq_5,Xi_eq_5,eye(M),gamma);
            mse_x_dsgs_5(node1,j) = norm(x-H*s_dsgs_5)^2;
            
            % Assuming BGS (bandlimited)- selecting K greatest coefficients
            [~, indexes] = sort(abs(V'*x),'descend');
            max_ind = sort(indexes(1:M));
            E_k = zeros(N,M);
            E_k_aux = eye(N,M);
            E_k(max_ind,:) = E_k_aux(1:M,:);
            Psi_i = Psi*diag(V(node1,:))*E_k;
            y_mi = Cm*Psi_i*V(:,max_ind)'*x;
            mse_bgs(node1,j) = norm(x-V(:,max_ind)*pinv(Cm*Psi_i)*y_mi)^2;
        end
    end

    g_prctile_mse(i,:,1) = prctile(mse_x_dsgs_1,percentile);
    g_prctile_mse(i,:,2) = prctile(mse_x_dsgs_2,percentile);
    g_prctile_mse(i,:,3) = prctile(mse_x_dsgs_5,percentile);
    g_prctile_mse(i,:,4) = prctile(mse_bgs,percentile);
end
time = toc;
disp(['Time: ' num2str(time/60) ' min'])

prctile_mse = zeros(length(MM),size(g_prctile_mse,3));
for k=1:size(g_prctile_mse,3)
    prctile_mse(:,k) = median(g_prctile_mse(:,:,k));
end

% Plot results
legend_txt = {'$$\hat{x}_{DSGS,1}$$', '$$\hat{x}_{DSGS,2}$$',...
              '$$\hat{x}_{DSGS,5}$$', '$$\hat{x}_{BGS}$$'};
figure();
semilogy(MM,prctile_mse(:,1),'-o','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(MM,prctile_mse(:,2),'-+','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(MM,prctile_mse(:,3),'-^','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(MM,prctile_mse(:,4),'-.*','LineWidth',1.5,'MarkerSize',8);hold off;
ylabel(['(b) MSE'],'FontSize', 20);
xlabel('Number of Observations','FontSize', 20);grid on;axis tight
legend(legend_txt, 'FontSize', 11,'Interpreter','latex');set(gca,'FontSize',14);
set(gcf, 'PaperPositionMode', 'auto')
