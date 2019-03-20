function [Cm_eq, Xi_eq] = eq_samp_matrices_k_nodes(N,M,nodes,Psi,V,h_freq)
% This function create the equivalent samplling matrices Cm and Xi for
% aggregation sampling using k nodes as the sampling nodes.

    % Divide the N and M nodes and observations on k equally spaced blocks
    k = length(nodes);
    Ns = divide_in_k_blocks(N,k);
    Ms = divide_in_k_blocks(M,k);

    Xi_eq = zeros(N);
    eyes = cell(k,1);
    cont = 1;
    for i=1:k
       eyes{i} = eye(Ms(i),Ns(i));
       Xi_i = Psi'*diag(V(nodes(i),:)'.*h_freq)*V';
       Xi_eq(cont:cont+Ns(i)-1,:) = Xi_i(1:Ns(i),:);
       cont = cont + Ns(i);
    end
    Cm_eq = blkdiag(eyes{:});
end

function Ns = divide_in_k_blocks(n,k)
% This function return a k length array with the sizes of the k blocks 
% resultng from dividing the block of size n as equal as posible

    Ns = zeros(k,1);
    for i=1:k
       if mod(n,k) == 0
           Ns(i) = n/k;
       else
           Ns(i) = ceil(n/k);
           n = n - 1;
       end
    end
end