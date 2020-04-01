function [nodes, fibers] = make_del(points_xyz)
% points_xyz -- N x 3 for N points
% x_col, y_col z_col -- column vectors
% ourmesh -- M x 4 node array for M tets corresponding to nodes array
% 
% last update -- mon aug 27 2012 -- mfh


x_col = points_xyz(:,1);
y_col = points_xyz(:,2);
z_col = points_xyz(:,3);

rand('state',sum(100*clock)); % scramble rand()

ourmesh = delaunay(x_col, y_col, z_col); % generate delaunay mesh
    
fibers = find_edges( ourmesh ); % find unique edges in mesh
    
nodes = points_xyz;

end
