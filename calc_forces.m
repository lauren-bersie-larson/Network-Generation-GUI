function [forces] = calc_forces(nodes, fibers, init_lens, modulii, fiber_area, fiber_B_vector)

% [fx] = calc_forces(n, f, i, m, a, b)
%
% calcs the acting fiber forces on nodes based on an exp fib exn -- in netmat
%
% nodes n -- 1 x 3N xyzxyz... nodal coordinates for N nodes
% fibers f -- 1 x 2N array of ababab... fiber nodes for N fibers
% init_lens i -- 1 x N initial lengths for N fibers
% modulii m -- 1 x N modulii for N fibers
% fiber_area a -- 1 x N cross-sectional areas for N fibers
% fiber_B_vector b -- 1 x N exp const exn b params for N fibers
%
% forces fx -- 1 x 3N array of xyzxyz... force components acting on N nodes
%
% fiber force is linearized above a critical value using an internal param
% this avoids the force becoming huge and generating a nan
% const exn for fiber forced used is F = E A (exp(B GS) - 1) / B
%
% depends on -- nothing
%
% last update -- sat july 28 2012 -- mfh


num_fibers = length(fibers) / 2;

nodes_x = nodes(1:3:end);
nodes_y = nodes(2:3:end);
nodes_z = nodes(3:3:end);

forces = zeros(1, length(nodes)); 

for n = 1 : num_fibers

    node_1_num = fibers(n*2-1);
    node_2_num = fibers(n*2-0);
    
    x_span = nodes_x( node_2_num ) - nodes_x( node_1_num );
    y_span = nodes_y( node_2_num ) - nodes_y( node_1_num );
    z_span = nodes_z( node_2_num ) - nodes_z( node_1_num ); 
    
    fiber_current_length = sqrt( x_span^2 + y_span^2 + z_span^2 );
   

    % EXP FIBER FORCE RELATION F = E A ( exp(B GS) - 1 ) / B
    
    lambda_limit = 4.0; % fib force linear above fiber green strain of 7.5
        
    lambda = fiber_current_length / init_lens(n);
    
    
    if lambda > lambda_limit

        % LINEAR FORCE ABOVE OUR THRESHOLD STRETCH
        
        gs = 0.5*(lambda_limit^2-1);
        
        force_exp = modulii(n) * fiber_area(n) * ...
                    (exp(fiber_B_vector(n) * gs)-1) / fiber_B_vector(n);
        
        slope_at_lambda_limit = modulii(n) * ...
                                fiber_area(n) * lambda_limit * ...
                                exp(fiber_B_vector(n) * gs);
        
        force = force_exp + slope_at_lambda_limit * (lambda - lambda_limit);
        
    else

        % REGULAR CONST EXN BELOW THRESHOLD
        
        gs = 0.5 * ( lambda^2 - 1 );
        
        force = modulii(n) * fiber_area(n) * ...
                (exp( fiber_B_vector(n) * gs) - 1) / fiber_B_vector(n);
        
    end
    
    cosine_A = x_span / fiber_current_length;
    cosine_B = y_span / fiber_current_length;
    cosine_C = z_span / fiber_current_length;

    node_1_force_x = force * cosine_A * (+1);
    node_1_force_y = force * cosine_B * (+1);
    node_1_force_z = force * cosine_C * (+1);

    node_2_force_x = force * cosine_A * (-1);
    node_2_force_y = force * cosine_B * (-1);
    node_2_force_z = force * cosine_C * (-1);

    forces( node_1_num*3-2 ) = forces(node_1_num*3-2) + node_1_force_x; 
    forces( node_1_num*3-1 ) = forces(node_1_num*3-1) + node_1_force_y; 
    forces( node_1_num*3-0 ) = forces(node_1_num*3-0) + node_1_force_z;

    forces( node_2_num*3-2 ) = forces(node_2_num*3-2) + node_2_force_x; 
    forces( node_2_num*3-1 ) = forces(node_2_num*3-1) + node_2_force_y; 
    forces( node_2_num*3-0 ) = forces(node_2_num*3-0) + node_2_force_z;

end


end
