function edges = find_edges(ourmesh)


% find edges in del mesh -- in netmat
%
% last update -- mon aug 27 2012 -- mfh


edges = []; % empty matrix to fill with polygon edges

poly_n = size( ourmesh, 1 ); % number of polygons

vertex_n = size( ourmesh, 2 ); % number of vertices per polygon



for n = 1 : poly_n

    nodes = ourmesh( n , : ); % nodes in polygon n
    
    edges = [ edges; nchoosek( nodes, 2 )]; % 2 nodes per edge
       
end


edges = sort( edges, 2 ); % sort nodes in each edge

edges = sortrows( edges ); % sort out edges by nodes just to look tidier

edges = unique( edges, 'rows' ); % our unique edges


end

