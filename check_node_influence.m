% Compare different nodes for sampling. Try to identify if selecting the
% node belonging to the same community as the data improves the resampling
addpath('./utils')

% Control variables
n_graphs = 500;
threshold = 0.01;
N = 50;
B = 5;
L = 6;
p = 0.4;
q = .1;
SS = 1:6;
M = 8;
p_n = 1e-5; %Noise power in natural units
gamma = 1e-1;
comm_nodes = B:B:N;

rec = zeros(10,length(SS));
tic
for k=1:length(SS)
    S = SS(k);
    fprintf('S: %d\n', S)

    mse_comm_known = zeros(n_graphs,length(comm_nodes));
    mse_comm_filt_known = zeros(n_graphs,length(comm_nodes));
    mse_comm_unknown = zeros(n_graphs,length(comm_nodes));
    mse_comm_filt_unknown = zeros(n_graphs,length(comm_nodes));
    mse_rand_known = zeros(n_graphs,length(comm_nodes));
    mse_rand_filt_known = zeros(n_graphs,length(comm_nodes));
    mse_rand_unkown = zeros(n_graphs,length(comm_nodes));
    mse_rand_filt_unkown = zeros(n_graphs,length(comm_nodes));
    mse_ss = zeros(n_graphs,length(comm_nodes));
    mse_ss_filt = zeros(n_graphs,length(comm_nodes));
    for i=1:n_graphs
        fprintf('\tGraph: %d\n', i)
        % Generate sparse signal s and filter coefficients h
        s = zeros(N,1);
        s(randperm(N/B,S)*B) = randn(S,1);
        Cs = zeros(S,N);
        non_zero_index = find(s~=0);
        for j=1:S
            Cs(j,non_zero_index(j)) = 1;
        end
        h = [randn(L,1); zeros(N-L,1)];
    
        % Generate graph
        [G, connected] = gen_est_conn_graph(N, B, p, q);
        if ~connected
            error('ERROR: graph is not connected')
        end
        Psi = fliplr(vander(G.e));
        
        % Sampling matrices and nodes list
        Cm = eye(M,N);
        Cms = rand_sampling_matrices(N,M,length(comm_nodes));
        nodes = 1:N;
        nodes = nodes(randperm(N,length(comm_nodes)));
        parfor j=1:length(comm_nodes)
            % Using AGSS with a node from the community
            comm_node = comm_nodes(j);
            Theta_i = Psi'*diag(G.U(comm_node,:))*G.U';
            Xi_i = Psi'*diag(G.U(comm_node,:)'.*(Psi*h))*G.U';
            wi = randn(N,1)*sqrt(p_n);
            y_mi = Cm*(Theta_i*s + wi);
            y_filt_mi = Cm*(Xi_i*s + wi);
            
            s_known = Cs'*pinv(Cm*Theta_i*Cs')*y_mi;
            s_filt_known = Cs'*pinv(Cm*Xi_i*Cs')*y_filt_mi;
            s_unk = estimate_noisy_s(y_mi,Cm,Theta_i,eye(M),gamma);
            s_filt_unk = estimate_noisy_s(y_filt_mi,Cm,Xi_i,eye(M),gamma);
            mse_comm_known(i,j) = norm(s-s_known)^2;
            mse_comm_filt_known(i,j) = norm(s-s_filt_known)^2;
            mse_comm_unknown(i,j) = norm(s-s_unk)^2;
            mse_comm_filt_unknown(i,j) = norm(s-s_filt_unk)^2;
            
            % Using AGSS with a random node
            rand_node = nodes(j);
            Theta_i = Psi'*diag(G.U(rand_node,:))*G.U';
            Xi_i = Psi'*diag(G.U(rand_node,:)'.*(Psi*h))*G.U';
            wi = randn(N,1)*sqrt(p_n);
            y_mi = Cm*(Theta_i*s + wi);
            y_filt_mi = Cm*(Xi_i*s + wi);
            
            s_known = Cs'*pinv(Cm*Theta_i*Cs')*y_mi;
            s_filt_known = Cs'*pinv(Cm*Xi_i*Cs')*y_filt_mi;
            s_unk = estimate_noisy_s(y_mi,Cm,Theta_i,eye(M),gamma);
            s_filt_unk = estimate_noisy_s(y_filt_mi,Cm,Xi_i,eye(M),gamma);     
            mse_rand_known(i,j) = norm(s-s_known)^2;
            mse_rand_filt_known(i,j) = norm(s-s_filt_known)^2;
            mse_rand_unkown(i,j) = norm(s-s_unk)^2;
            mse_rand_filt_unkown(i,j) = norm(s-s_filt_unk)^2;
            
            % Using SS - At least one node is from the same community as s
            Cm_ss = Cms(:,:,j);
            H = G.U*diag(Psi*h)*G.U';
            x = H*s;
            s_ss = estimate_noisy_s(Cm_ss*s,Cm_ss,eye(N),eye(M),gamma);
            s_ss_filt = estimate_noisy_s(Cm_ss*x,Cm_ss,H,eye(M),gamma);
            mse_ss(i,j) = norm(s-s_ss);
            mse_ss_filt(i,j) = norm(s-s_ss_filt)^2;
        end
    end

    rec(1,k) = median(sum(mse_comm_known < threshold)./n_graphs);
    rec(2,k) = median(sum(mse_comm_filt_known < threshold)./n_graphs);
    rec(3,k) = median(sum(mse_comm_unknown < threshold)./n_graphs);
    rec(4,k) = median(sum(mse_comm_filt_unknown < threshold)./n_graphs);
    rec(5,k) = median(sum(mse_rand_known < threshold)./n_graphs);
    rec(6,k) = median(sum(mse_rand_filt_known < threshold)./n_graphs);
    rec(7,k) = median(sum(mse_rand_unkown < threshold)./n_graphs);
    rec(8,k) = median(sum(mse_rand_filt_unkown < threshold)./n_graphs);
    rec(9,k) = median(sum(mse_ss < threshold)./n_graphs);
    rec(10,k) = median(sum(mse_ss_filt < threshold)./n_graphs);
end
time = toc;
disp(['Time: ' num2str(time/60) ' min'])

% Plot results
legend_txt = {'$$\hat{s}_{AGSS,Comm,Known}$$','$$\hat{x}_{AGSS,Comm,Known}$$',...
              '$$\hat{s}_{AGSS,Comm,Unk}$$','$$\hat{x}_{AGSS,Comm,Unk}$$',...
              '$$\hat{s}_{AGSS,Rand,Known}$$', '$$\hat{x}_{AGSS,Rand,Known}$$',...
              '$$\hat{s}_{AGSS,Rand,Unk}$$', '$$\hat{x}_{AGSS,Rand,Unk}$$',...
              '$$\hat{s}_{SS,Unk}$$', '$$\hat{x}_{SS,Unk}$$'};          
figure()
plot(SS,rec(1,:),'-+','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(2,:),'-*','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(3,:),'-^','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(4,:),'-o','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(5,:),'-.+','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(6,:),'-.*','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(7,:),'-.^','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(8,:),'-.o','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(9,:),'-s','LineWidth',1.5,'MarkerSize',8);hold on
plot(SS,rec(10,:),'-.s','LineWidth',1.5,'MarkerSize',8);
xlabel('Number of Seeding Nodes','FontSize', 20)
ylabel('(a) Recovery Rate','FontSize', 20)
legend(legend_txt, 'FontSize', 14,'Interpreter','latex')
grid on; set(gcf, 'PaperPositionMode', 'auto'); set(gca,'FontSize',14);
