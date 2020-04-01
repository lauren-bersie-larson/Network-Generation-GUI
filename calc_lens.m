function [lens] = calc_lens(nodes, fibers)
   global dim

    % calculate length of network fibers for 1D inputs -- in netmat
    % nodes -- 1 x 3N vect of xyz coords for N nodes
    % fibers -- 1 x 2N vect of start-end node nums for N fibers
    % lens -- 1 x N vect of N fiber lengths
    
    % last update -- mon jul 23 2012 -- mfh
    % Updated 6-26-17 LMB to include dim switch- for 2D vs 3D networks

    switch dim
        case 2
            num_fibers = length(fibers) / 2;
            lens = zeros(1, num_fibers);

            nodes_x = nodes(1:2:end);
            nodes_y = nodes(2:2:end);

            for n = 1 : num_fibers
                node_1_num = fibers(n*2-1);
                node_2_num = fibers(n*2-0);
    
                node_1_vect=[nodes_x(node_1_num) nodes_y(node_1_num)];
                node_2_vect=[nodes_x(node_2_num) nodes_y(node_2_num)];
    
                lens(n) = norm(node_1_vect - node_2_vect);
            end
        case 3
            num_fibers = length(fibers) / 2;
            lens = zeros(1, num_fibers);

            nodes_x = nodes(1:3:end);
            nodes_y = nodes(2:3:end);
            nodes_z = nodes(3:3:end);

            for n = 1 : num_fibers
                node_1_num = fibers(n*2-1);
                node_2_num = fibers(n*2-0);
    
                node_1_vect=[nodes_x(node_1_num) nodes_y(node_1_num) nodes_z(node_1_num)];
                node_2_vect=[nodes_x(node_2_num) nodes_y(node_2_num) nodes_z(node_2_num)];
    
                lens(n) = norm(node_1_vect - node_2_vect);
            end
    end
end
