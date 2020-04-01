function [nodes_connectivity] = calc_degree(nodes, fibers)


% [nodes_connectivity] = calc_degree(nodes, fibers)
%
% calculate the connectivity of nodes in a network
%
% in:
% nodes                 N x 3 for N nodes
% fibers                M x 2 for M fibers
%
% out:
% nodes_connectivity    N x 1 for N nodal degrees


num_fibers = size(fibers, 1); % number of rows = number of fibers

num_nodes = size(nodes, 1); % number of rows = number of nodes


nodes_connectivity = zeros(num_nodes, 1); 

% we will return a column vector -- row number = node number

for n = 1 : num_fibers
    
    node1 = fibers(n,1);
    node2 = fibers(n,2);
    
    nodes_connectivity( node1 , 1 ) = nodes_connectivity( node1 , 1 ) + 1;
    nodes_connectivity( node2 , 1 ) = nodes_connectivity( node2 , 1 ) + 1;
    
end

end
