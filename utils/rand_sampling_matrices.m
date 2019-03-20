function [Cms, nodes] = rand_sampling_matrices(N,M, n_matrices)
% n_matrices is assumed to be less or equal to N
% For the fisrt n_matrices nodes, each node is forced to be a sampling node
% at least in one of the combiations

nodes = zeros(M,n_matrices);
Cms = zeros(M,N,n_matrices);

if n_matrices > (nchoosek(N,M)-1)
    disp('ERROR: there are not enough combinations')
    return
end

if n_matrices > N
    disp('ERROR: n_matrices is supposed to be <= than N')
    return
end

for i=1:n_matrices    
    nodes_aux = 1:N;
    nodes_aux(i) = [];
    nodes(:,i) = sort([i nodes_aux(randperm(N-1,M-1))]);
    while sum(ismember(nodes(:,1:i-1)',nodes(:,i)','rows')) ~= 0
        nodes(:,i) = sort([i nodes_aux(randperm(N-1,M-1))]);
    end
    for j=1:M
        Cms(j,nodes(j,i),i) = 1;
    end
end    