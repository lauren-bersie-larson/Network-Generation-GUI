function [segments] = find_connected(nodes, fibers)


% [segments] = find_connected(nodes, fibers)
%
% finds connected fiber segment groups in networks -- in netmat
%
% segments -- structure [segments(J).nodes] of connected segments
%
% nodes -- N x 3 nodal xyz coordinate rows for N nodes
% fibers -- N x 2 start-end nodes for N fibers
%
% depends on -- calc_adjacency(), dmperm()
% 
% last update -- fri aug 17 2012 -- mfh


num_fibers  = size(fibers,1);   % rows = num fibers

num_nodes   = size(nodes,1);    % rows = num nodes


% calculate the adjacency matrix

[adjacency_matrix] = calc_adjacency(fibers, nodes);



% next lump nodes into connected sets
%
% We'll use the Dulmage-Mendelsohn decomposition of a matrix to do this.
%
% This method has been taken directly from:
%
% blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/


% We'll start by filling the diagonal of the adjacency matrix with 1's

adjacency_matrix( 1 : num_nodes+1 : end ) = 1; % fill diagonal with 1's



% Next we'll call DMPERM() the Dulmage-Mendelsohn decomposition

[p,q,r,s] = dmperm(adjacency_matrix);


% [P] contains the full node set
% [R] contains indices which point to connected segments of the network



% We'll assume we have at least ONE connected segment

num_connected_segments = length(r) - 1;



for n = 1 : num_connected_segments 
    
    start_index = r(n);         % find start index from r
    
    end_index = r(n+1) - 1;     % find end index from r
   
    segments(n).nodes = p( start_index : end_index );
      
end


end
