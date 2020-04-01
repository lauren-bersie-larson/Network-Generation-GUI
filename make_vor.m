function [nodes, fibers] = make_vor(points_xyz)


% make_vor(points_xyz)
%
% generate a voronoi network from a set of seed points -- in netmat
%
% depends on -- voronoin(), convhulln(), TriRep()
%
% contains -- filter_cells(), renum_nodes()
%
% updated -- sat aug 18 2012 -- mfh


% turn seed points into row vectors of x y z values

x = points_xyz(:,1)';
y = points_xyz(:,2)';
z = points_xyz(:,3)';


% apply voronoi tesselation to points

[nodes, cells] = voronoin( [x(:),y(:),z(:)] ); 



% remove all cells connected to node 1 which is at <inf>

cells = filter_cells( cells , nodes);



% generate point cloud for each voronoi cell

sharp_edges = [];


for n = 1 : length( cells )
    
   unique_vertices_in_bubble = cells{n};
   
   node_xyz = nodes( unique_vertices_in_bubble, : );
  
 
   % CONVHULLN returns the triangular facets of the cube faces
   
   node_vertices = convhulln(node_xyz);

   
   % TRIREP is used to find the FEATUREDGES of the cube faces
   
   node_tri = TriRep(node_vertices, node_xyz);
   
   
   local_edges = featureEdges(node_tri, pi/1000000); % <== threshold for surf
   
   
   
   % CONVERT LOCAL TO GLOBAL INDEXING
   
   for i = 1 : size( local_edges, 1 )
       
      global_edges(i,1) =  unique_vertices_in_bubble( local_edges(i,1) );
      
      global_edges(i,2) =  unique_vertices_in_bubble( local_edges(i,2) );
       
   end
   
   
   sharp_edges = [ sharp_edges;
                   global_edges ]; 

    
end



% GATHER UNIQUE FIBERS

fibers = unique(sharp_edges, 'rows');



% RE-NUMBER NETWORK NODES

[nodes, fibers] = renum_nodes(nodes, fibers);



end



function [new_cells] = filter_cells(cells, nodes)

% nodes -- m-by-3 array of xyz coordinates for m nodes indexed 1-m
% fibers -- m-by-2 array of start-end node indices for 1-m fibers
%
% new_cells -- all cells except ones with node 1 numbered 1...N cells
%
% updated -- sun aug 19 2012 -- mfh

new_cell_count = 1;

for n = 1 : length( cells )  
   
    if cells{n}(1) ~= 1; % Avoid the cells that have the infinite node number 1
        
        new_cells{ new_cell_count } = cells{ n };
        
        new_cell_count = new_cell_count + 1;
       
    end
    
end

end



function [new_nodes, new_fibers] = renum_nodes(nodes, fibers)


% Function writes out a network where some instances in NODES do not
% appear in the FIBERS array -- such as in the case of a Voronoi network
% where we skip cells linked to an infinite node -- renumbering the nodes
% that actually do appear in the network.
%
% INPUT
% =====
% NODES - m-by-3 array of xyz coordinates for m nodes indexed 1-m
% FIBERS - m-by-2 array of start-end node indices for 1-m fibers
% FILENAME - string containing filename to write to
%
% OUTPUT
% ======
% NODES - resized to remove unused nodes
% FIBERS - nodes re-numbered sequentially
%
% HISTORY
% =======
% Apr 2011 -- Modified from write_network() -- MFH



total_fibers = size( fibers, 1 ); % num rows=num fibers regardless of indexing



% Collect old node indices
% ------------------------

old_node_indices = unique( fibers(1:end) ); % 1 x N array -- sorted

total_nodes = length( old_node_indices ); % number of nodes in play




% Renumber FIBERS node indices and collect new NODES array
% --------------------------------------------------------

for n = 1 : total_fibers
    
   old_start_node       = fibers(n,1);
   old_end_node         = fibers(n,2);
   
   new_start_index      = find( old_node_indices == old_start_node );
   new_end_index        = find( old_node_indices == old_end_node );
   
   fibers(n,1)          = new_start_index;
   fibers(n,2)          = new_end_index;
    
end


new_nodes = nodes(old_node_indices, :); % return new nodes

new_fibers = fibers; % just to be clear we updated fibers to return

end




