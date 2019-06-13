function [X,colors,known_nodes]=etex_measures2array(measures_data,debug)
% Constants
GREEN = [0 1 0];
RED = [1 0 0];
BLUE = [0 0 1];

% Data Flags
NOT_SAMP = -0.99;
NOT_QUANT = -0.88; 

known_nodes = cell(length(measures_data)-2,1);

colors = repmat(BLUE,[length(measures_data{1}),1,30]);
X = zeros(length(measures_data{4}),length(measures_data)-2);
for i=3:length(measures_data)
    X(:,i-2) = measures_data{i};
    ind_measured = X(:,i-2)>0;
    ind_not_samp = X(:,i-2)==NOT_SAMP;
    ind_not_quant = X(:,i-2)==NOT_QUANT;
    
    if debug
        disp(['Summary X[',num2str(i-2), ']:'])
        tr_measured = sum(ind_measured)/length(X(:,i-2))*100;
        not_samp = sum(ind_not_samp)/length(X(:,i-2))*100;
        not_quant = sum(ind_not_quant)/length(X(:,i-2))*100;
        disp(['   Tracer meassured: ',num2str(tr_measured),'%   Not samp: ',...
                num2str(not_samp),'%   Not quant: ', num2str(not_quant) '%'])
    end
 
    colors(ind_measured,:,i-2)=repmat(GREEN,[sum(ind_measured) 1 1]);
    colors(ind_not_samp|ind_not_quant,:,i-2)=repmat(RED,...
                                    [sum(ind_not_samp|ind_not_quant) 1 1]);
    known_nodes{i-2} = find(X(:,i-2)>=0);
end

X(X<0) = 0;