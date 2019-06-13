function [A,X] = remove_nodes(A,X,disc_nodes)
A(disc_nodes,:)=[];
A(:,disc_nodes)=[];

for i=1:length(X)
    X{i}(disc_nodes)=[];
end
