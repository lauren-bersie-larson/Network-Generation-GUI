function [new_nodes, new_fibers] = get_giant(nodes, fibers)


% [new_nodes, new_fibers] = get_giant(nodes, fibers)
%
% returns giant network within a fiber network -- in netmat
%
% nodes -- N x 3 input nodes
% fibers -- N x 2 input fibers
%
% new_nodes -- N x 3 giant net nodes
% new_fibers -- N x 2 giant net fibers
%
% depends on -- find_connected(), cellfun(), max()
%
% last update -- tue aug 21 2012 -- mfh


% FIND CONNECTED REGIONS AND NODES TO KEEP

[segments] = find_connected(nodes, fibers);

cellfun(@length, {segments.nodes});

% find largest connected segment -- its index is the one to target
[mag,ind] = max(cellfun(@length, {segments.nodes}));

node_nums_to_keep = [segments(ind).nodes]; % pulls out largest connectome


% ADD EXTRA COLS FOR NODE AND FIBER NUMS

nodes = [ [1:size(nodes,1)]' nodes];

fibers = [ [1:size(fibers,1)]' fibers];



% FIND FIBERS TO KEEP

[not_unique_fiber_rows_to_keep, not_unique_node_cols] = ...
            find(ismember(fibers(:, 2:3), node_nums_to_keep));

fiber_rows_to_keep = unique(not_unique_fiber_rows_to_keep);

fibers_to_keep = fibers(fiber_rows_to_keep, :);

nodes_to_keep = nodes(node_nums_to_keep, :);


% CLEAN UP FIBERS AND NODES TO REMOVE NUMBERING GAPS

[final_nodes, final_fibers] = giant_net_cleanup(nodes_to_keep, fibers_to_keep);


new_nodes = final_nodes(:, 2:4);

new_fibers = final_fibers(:, 2:3);


end





function [final_nodes, final_fibers] = giant_net_cleanup(new_node_stack, new_fiber_stack)


% NET_CLEANUP
%
% new_node_stack   4 x N col one is node number col two-four node xyz
% new_fiber_stack  3 x M col one is fiber number col two-three start-end nodes
%
% final_nodes      4 x I col one is node number col two-four node xyz
% final_fibers     3 x J col one is fiber number col two-three start-end nodes


% net_cleanup() should be able to handle gaps in node numbering and 
% fiber numbering so that we can remove fibers and nodes


% SORT OUT AND RENUMBER NODES
% ---------------------------

new_node_stack = sortrows( unique( new_node_stack, 'rows')); 

% UNIQUE() should also sort in
% order via the first column of
% values (ie node number) but
% we'll re-sort just in case

                                                             
number_of_nodes = size(new_node_stack, 1); % number of rows = number nodes



new_node_list = 1 : number_of_nodes; 

old_node_list = new_node_stack( : , 1 );


new_node_stack( : , 1 ) = 1 : number_of_nodes; % renumber nodes



% RENUMBER FIBER START-END NODES
% ------------------------------


number_of_fibers = size( new_fiber_stack, 1 ); % count rows = num fibers

for n = 1:number_of_fibers
    
   old_node_a = new_fiber_stack( n , 2 );
   old_node_b = new_fiber_stack( n , 3 );
   
   old_node_index_a = find( old_node_list == old_node_a );
   old_node_index_b = find( old_node_list == old_node_b );
   
   new_node_a = new_node_list( old_node_index_a );
   new_node_b = new_node_list( old_node_index_b );
   
   new_fiber_stack( n , 2 ) = new_node_a; 
   new_fiber_stack( n , 3 ) = new_node_b;
    
end

% renumber fibers

new_fiber_stack(:,1) = 1 : number_of_fibers; % first col = fiber number


% output resultant fibers

final_nodes     =   new_node_stack;   % nodes cleaned up
final_fibers    =   new_fiber_stack;  % fibers cleaned up


end
