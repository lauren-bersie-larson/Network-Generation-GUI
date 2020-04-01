function [adjacency_matrix] = calc_adjacency(fibers, nodes)


% [adjacency_matrix] = calc_adjacency(fibers, nodes)
%
% returns the adjacency matrix (zero diagonal) for a network -- in netmat
%
% nodes -- N x 3 nodal xyz coordinate rows for N nodes
% fibers -- N x 2 start-end nodes for N fibers
%
% updated -- fri aug 17 2012 -- mfh


num_fibers  = size(fibers,1);   % rows = num fibers
num_nodes   = size(nodes,1);    % rows = num nodes

adjacency_matrix = zeros(num_nodes, num_nodes);

% populate adjacency matrix

for n = 1 : num_fibers
    
   node_A = fibers(n,1);
   node_B = fibers(n,2);
   
   adjacency_matrix(node_A, node_B) = 1;
   adjacency_matrix(node_B, node_A) = 1; 
    
end

end
