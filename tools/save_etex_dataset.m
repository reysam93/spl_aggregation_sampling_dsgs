% Constants
STATIONS_PATH = '..\dataset_additional\stationlist.950130';
MEASURES_PATH = '..\dataset_additional\pmcp.datv1.2';

% Variables
max_km = 200;
debug_dataset = true;
save_data = false;
gso = 'A';
geo_graph = false;

% Read original data
fid = fopen(STATIONS_PATH);
stations_data = textscan(fid,'%s%s%f%f%d%d%*[^\n]',168,'Delimiter',',','HeaderLines',5);
fclose(fid);
fid = fopen(MEASURES_PATH);
measures_data = textscan(fid,'%d%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'...
                        ,'HeaderLines',2);
fclose(fid);

% Construc graph, signals and node labels
if geo_graph
    DISCON_NODES = [14 15 134 136 137];
    [A,D] = adajacency_from_stations(stations_data, max_km);
    [A,measures_data] = remove_nodes(A,measures_data,DISCON_NODES);
else
   load('..\dataset_additional\etex_data_from_paper.mat');
   A=double(full(G.A));
end
[X, node_colors,known_nodes] = etex_measures2array(measures_data,debug_dataset);

% Plot Graph
if debug_dataset
    figure();spy(A)
    G = graph(A);
    figure();
    subplot(1,3,1)
    plot(G,'MarkerSize',4,'NodeColor',node_colors(:,:,1));
    subplot(1,3,2)
    plot(G,'MarkerSize',4,'NodeColor',node_colors(:,:,15));
    subplot(1,3,3)
    plot(G,'MarkerSize',4,'NodeColor',node_colors(:,:,25));
end


%% Estimate Filter
s = X(:,1);
n_graphs = size(X,2)-1;
hh_freq = cell(n_graphs,1);
vv = cell(n_graphs,1);
ee = cell(n_graphs,1);
err_no_skip = zeros(n_graphs,1);
mean_err = zeros(n_graphs,1);
mean_sparsity = zeros(n_graphs,1);
skipped = 0;
hh = cell(n_graphs,1);
c_knowns = cell(n_graphs,1);
for i=1:n_graphs
    x = X(known_nodes{i+1},i+1);
    hh_freq{i} = [];
    
    % if x is 0 --> skip!
    if sum(x)==0
        skipped = skipped + 1;
        continue
    end
    
    % Eigendecomposition of the GSO
    if strcmp(gso,'L')
        G = gsp_graph(A);
        G = gsp_compute_fourier_basis(G);
        e = G.e;
        V = G.U;
        GSO = G.L;
    elseif strcmp(gso,'A')
        [V,D] = eig(A);
        e = diag(D);
        GSO = A;
    end
    
    Cm = zeros(length(x),length(A));
    for j=1:length(x)
        Cm(j,known_nodes{i+1}(j))=1;
    end
    c_knowns{i} = Cm;
    
    h_freq = estimate_h_freq_known_s(x,s,V,Cm);
    if isempty(h_freq)
        skipped = skipped + 1;
        continue
    end
    hh_freq{i} = h_freq;
    
    vv{i} = V;
    ee{i} = e;
    err_no_skip(i) = norm(Cm*V*diag(hh_freq{i})*V'*s-x);
    disp(['Graph ' num2str(i+1) '   N: ' num2str(length(x)) '   Sparsity: '... 
        num2str(sum(abs(s)>0)) '   MSE: ' num2str(err_no_skip(i))])
end

% Remove skipped graphs from data
SS = cell(n_graphs-skipped,1);
XX = cell(n_graphs-skipped,1);
HH_Freq = cell(n_graphs-skipped,1);
VV = cell(n_graphs-skipped,1);
EE = cell(n_graphs-skipped,1);
C_KNOWNS = cell(n_graphs-skipped,1);
Time_Index = zeros(n_graphs-skipped,1); 
err_skip = zeros(n_graphs-skipped,1);
max_h = zeros(n_graphs-skipped,1);
cont = 1;
for j=1:(n_graphs)
    if isempty(hh_freq{j})
        continue
    end
    XX{cont} = X(known_nodes{j+1},j+1);
    C_KNOWNS{cont} =  c_knowns{j};
    SS{cont} = s;
    HH_Freq{cont} = hh_freq{j};
    VV{cont} = vv{j};
    EE{cont} = ee{j};
    Time_Index(cont) = j;
    max_h(cont) = max(hh_freq{j});
    err_skip(cont) = err_no_skip(j);
    cont = cont + 1;
end

disp(['Mean MSE:        ' num2str(mean(err_skip))])
disp(['Selected graphs: ' num2str(length(XX))])
disp(['Mean Max h       ' num2str(mean(max_h))])

if save_data
     save(['.\ETEX_datasets\LOGIC_GRAPH_UNK_etex1_eig' gso '_e-4'],...
            'EE', 'VV', 'XX', 'SS', 'HH_Freq', 'C_KNOWNS', 'Time_Index')
end


