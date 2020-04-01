function [] = put_net(nodes, fibers, fnm, fib_reg)


% [] = put_net(nodes, fibers, filename)
%
% writes out a network to a single text file -- in netmat
%
% filename -- string with path to file
% nodes -- N x 3 nodal xyz coordinate rows for N nodes
% fibers -- N x 2 start-end nodes for N fibers
% Edit: Stores file in Textfiles/ now
% last update -- Thu 13 Feb 2014 -- RD


total_nodes         = size(nodes, 1);     % count rows = num nodes
total_fibers        = size(fibers, 1);    % count rows = num fibers

tot_deg_freedom     = 3 * total_nodes;      % x y z values for nodes

filename = fnm;%fullfile('Textfiles', fnm);
fileid = fopen(filename , 'w'); % writes over any existing file

fprintf(fileid, '%i %i %i\n', total_nodes, tot_deg_freedom, total_fibers);

for n = 1 : total_fibers

    fprintf(fileid,'%i %i %i %f %f %f %f %f %f %i\n', ...  % file line filters
        n, ...                                          % fiber num
        fibers(n,1), ...                                % node 1 num
        fibers(n,2), ...                                % node 2 num
        nodes( fibers(n,1), 1 ), ...                    % node 1 x coord
        nodes( fibers(n,1), 2 ), ...                    % node 1 y coord
        nodes( fibers(n,1), 3 ), ...                    % node 1 z coord
        nodes( fibers(n,2), 1 ), ...                    % node 2 x coord
        nodes( fibers(n,2), 2 ), ...                    % node 2 y coord
        nodes( fibers(n,2), 3 ), ...                   % node 2 z coord
        fib_reg(n,1));                                   % Fiber type

end

fclose(fileid);

end

