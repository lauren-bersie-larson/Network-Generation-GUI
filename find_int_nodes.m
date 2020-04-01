function [int_node_nums] = find_int_nodes(nodes,boundaries)


% Out
% ===
% int_node_nums -- 1 x N vector of N interior node numbers


% SETUP VARIABLES

xmin = boundaries(1);
xmax = boundaries(2);
ymin = boundaries(3);
ymax = boundaries(4);
zmin = boundaries(5);
zmax = boundaries(6);

num_nodes = length(nodes) / 3;

nodes_x = nodes(1:3:end);
nodes_y = nodes(2:3:end);
nodes_z = nodes(3:3:end);

int_node_nums = []; % will populate with bnd node indices



% FIND BOUNDARY NODES GENERICALLY


for n = 1 : num_nodes
    
   if( nodes_x(n)~=xmin && nodes_x(n)~=xmax && nodes_y(n)~=ymin && nodes_y(n) ~= ymax && nodes_z(n)~=zmin && nodes_z(n) ~= zmax )
       
       int_node_nums = cat(2, int_node_nums, n); %[int_node_nums n];
       
   end    
    
end




end




