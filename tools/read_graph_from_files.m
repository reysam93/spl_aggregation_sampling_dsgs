function [cell_A, cell_X, g_labels] = read_graph_from_files(files, n_graphs)
% Read the graph data obtained from https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
% and return two cell arrays with the adyacency matrix and the node's attributes
% Arguments:
%   - files: a struct containing the file's path. Must have the fields 'A',
%           'graph_indicator' and 'node_attribute'
%   - n_graphs: the number of graph desired to read. If omitted, all graphs will be read 

cell_A = {};
cell_X = {};
if ~isfield(files,'A') || ~isfield(files,'graph_indicator')...
        || ~isfield(files,'node_attribute') || ~isfield(files,'graph_label')
    disp('Missing struct fields')
    return 
end

AA = dlmread(files.A);
g_indicator = dlmread(files.graph_indicator);
nodes_attr = dlmread(files.node_attribute);
g_labels = dlmread(files.graph_label);

if nargin < 2 || n_graphs > g_indicator(end)
   n_graphs = g_indicator(end);
end

cell_A = cell(n_graphs,1);
cell_X = cell(n_graphs,1);
first_node = 1;
for i=1:n_graphs
    n_nodes = sum(g_indicator==i);
    last_node = first_node+n_nodes-1;
    cell_X{i} = nodes_attr(first_node:last_node,:);
    
    A = zeros(n_nodes);
    idx = AA(:,1)>=first_node & AA(:,1)<=last_node;
    row_idx = AA(idx,1)-(first_node-1);
    col_idx = AA(idx,2)-(first_node-1);
    A(sub2ind(size(A),row_idx,col_idx)) = 1;
    
    assert(sum(diag(A))==0,'ERR: diagonal of A is not 0')
    assert(sum(sum(A-A'))==0,'ERR: A is not symmetric')
    
    cell_A{i} = A;
    first_node = last_node + 1;
end
