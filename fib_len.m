function [lens] = fib_len(nodes, fibers)

% [lens] = fib_len(nodes, fibers)
%
% calculate the length of fibers for 2D inputs -- in netmat
%
% nodes -- N x 3 nodal xyz coordinate rows for N nodes
% fibers -- N x 2 start-end nodes for N fibers
% lens -- N x 1 lengths for N fibers
%
% last update -- jul 2012 -- mfh

num_fibers = size(fibers, 1); % number of rows = number of fibers

lens = zeros(num_fibers, 1); % we will return a column vector

for n = 1 : num_fibers
    
    node1 = fibers(n,1); % node 1 num
    node2 = fibers(n,2); % node 2 num
   
    % fiber length = norm([x1 y1 z1] - [x2 y2 z2])
 
    lens(n) = norm( [nodes(node1,:)] - [nodes(node2,:)] );
    
end

end
