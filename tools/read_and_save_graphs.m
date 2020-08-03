% Read all the graphs from the text files and save all of them and the
% first 100 graphs
save_files = false;

files = struct();
files.A = '..\dataset_additional\PROTEINS_full_A.txt';
files.graph_indicator = '..\dataset_additional\PROTEINS_full_graph_indicator.txt';
files.graph_label = '..\dataset_additional\PROTEINS_full_graph_labels.txt';
files.node_attribute = '..\dataset_additional\PROTEINS_full_node_attributes.txt';
n_graphs = 100;

[cell_A, cell_X, g_labels] = read_graph_from_files(files); %whole dataset
if save_files
    save('..\dataset_additional\graphs','cell_A', 'cell_X','g_labels')
end

cell_A = cell_A(1:n_graphs);
cell_X = cell_X(1:n_graphs);

if save_files
    save(['..\dataset_additional\read_' num2str(n_graphs) '_graphs'],'cell_A',...
        'cell_X','g_labels')
end