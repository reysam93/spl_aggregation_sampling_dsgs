function [G, connected] = gen_est_conn_graph(N, B, p, q, directed)
% Generate a conencted graph using the estochastic block model. It
% always distrubute the number of nodes N homogenuosly among all
% communities B.
% Arguments:
%   - N: number of nodes
%   - B: number of communities
%   - p: probability of creating a link between nodes of the same cluster
%   - q: probability of creating a link between nodes of different clusters
%   - directed: whether the graph is directed or undirected
if ~exist('directed', 'var')
    directed = false;
end

max_retries = 20;
params.p = p;
params.q = q;
params.z = ones(1,N);

% Grapsh are only generated as directed graph with force_full=true
params.directed = directed;
params.force_full = directed;

% params.z: vector indicating for each node to which cluster belongs
% generating all the communities with the same number of nodes
for i=1:N
    params.z(i) = params.z(i)*(mod(i,B)+1);
end

G = gsp_stochastic_block_graph(N, B, params);
G = gsp_compute_fourier_basis(G);
retries = 0;
connected = false;
while (sum(G.e > 1e-4) < N-1)
    G = gsp_stochastic_block_graph(N, B, params);
    G = gsp_compute_fourier_basis(G);
    retries = retries + 1;
    if retries >= max_retries
        disp('Could not create a connected graph')
        return
    end
end
connected = true;
