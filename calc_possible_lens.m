function [A] = calc_possible_lens(nodes)
% calc all possible fiber lens from a set of nodes
% nodes         N x 3 for N nodes
% last update wed dec 19 2012

n_nodes = size(nodes, 1); % rows

A = zeros(n_nodes, n_nodes);

for n = 1 : n_nodes % start nodes
    for m = (n+1) : n_nodes % end nodes
        % fiber len from nth to mth node
        A(n, m) = norm( [nodes(n, :) - nodes(m, :)] );        
        A(m, n) = A(n, m); % symmetric 
    end
end
end