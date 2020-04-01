%%% Script to create networks with a certain volume fraction 
%%% Written by Lauren Bersie, 7-19-17
%%% Updated 3-1-2020 to have equations for quickly finding seed point numbers, 
%%% GUI control, and network rotations

% DEPENDENCIES: make_aligned_dels_LB, calc_lens, plot_net
function create_volfract(net_name, net_type, num_nets, rve_len, target_vol_fract, align_dir, net_align, net_rot, rot_axis, save_path, fiber_rad)
global prog_box points_seed dim N % N is the number of networks to create 

% Set parameters for network generation
dim = 3; % 3D networks
boundaries = [-0.5 0.5 -0.5 0.5 -0.5 0.5]; % RVE boundaries in computational units
vf_tol = 1e-3; % Tolerance for final volume fraction. This ~ volume of one fiber
x = rve_len * 1e-6; % Dimensionalize to be in m- currently um (Change later)
fiber_radius = fiber_rad * 1e-9; % Dimensionalize to be in m- currently nm
x_incr = x/10e-6; % Equations for seed points based on dim of 20e-6. Will need to scale accordingly
vf_incr = target_vol_fract/0.04; %Equations for seed points based on vf of 0.04. Will need to scale accordingly

% Calculate network stretch needed for desired alignment
if strcmp(align_dir, 'X')
        lambdax = 0.355*exp(3*net_align);
        lambday = 1;
        lambdaz = 1;

elseif strcmp(align_dir, 'Y')
        lambdax = 1;
        lambday = 0.355*exp(3*net_align);
        lambdaz = 1;

elseif strcmp(align_dir, 'Z')
        lambdax = 1;
        lambday = 1;
        lambdaz = 0.355*exp(3*net_align);
end
    
if strcmp(net_type, 'Delaunay')  
    % Calculate number of seed points needed for given volume fraction
    if strcmp(align_dir, 'X')
            points_seed = (-273*(lambday) + 2940) * vf_incr * x_incr;

    elseif strcmp(align_dir, 'Y')
            points_seed = (-273*(lambday) + 2940) * vf_incr * x_incr;

    elseif strcmp(align_dir, 'Z')
            points_seed = (-273*(lambdaz) + 2940) * vf_incr * x_incr;
    end
    
elseif strcmp(net_type, 'Voronoi')
    % Calculate number of seed points needed for given volume fraction
    if strcmp(align_dir, 'X')
            %points_seed = (-760*(lambdax) + 7080) * vf_incr * x_incr;
            points_seed = (5500) * vf_incr * x_incr;
    elseif strcmp(align_dir, 'Y')
            points_seed = (-760*(lambday) + 7080) * vf_incr * x_incr;

    elseif strcmp(align_dir, 'Z')
            points_seed = (-760*(lambdaz) + 7080) * vf_incr * x_incr;
    end
    
end

% vf_incr = 1; % Change this percentage increase to change amount of volume fraction in networks
% target_vol_fract = vf_incr * 0.04; 
% vf_tol = 1e-3; % Tolerance for meeting target vol fraction (around the vol of one fiber) CHANGE
% major_radius = 100e-9; % Radius of major chain fibers, in m 
% RVE_len = 20e-6; % RVE side length, in m   
% x_incr = 1; % % How much you want to scale the RVE side length by
% x = RVE_len * x_incr;

% Make networks
length_vect = zeros(1,7);
orientation = zeros(3,3,7);
%     if vf_incr ~= 1
    points_seed = round(points_seed); % Adjust initial guess to make faster
%     end
int_nodes = [];
num_int_nodes = [];

for j = 1:num_nets % Create j networks
    % Update progress bar 
    prog_box.Value = j/num_nets; 
    message = sprintf('Creating network %i of %i. . .', j, num_nets);
    prog_box.Message = message;

    N = j;
    
    % Create network based on initial # of seed points     
    if strcmp(net_type, 'Delaunay') 
        [nodes_old, nodes, fibers_old, fibers] = make_aligned_dels_LB();
    elseif strcmp(net_type, 'Voronoi')
        [nodes_old, nodes, fibers_old, fibers] = make_aligned_vors(); 
    end

    % Determine the volume fraction of the network and compare to target
    nodes1D = conv_2D_2_lin(nodes); % Reshape into 1D array
    fibers1D = conv_2D_2_lin(fibers); % Reshape into 1D array
    [lens] = calc_lens(nodes1D, fibers1D); % Returns lengths of every fiber
    total_len = sum(lens);
    % Calculate volume fraction    
    total_vf = (total_len*pi*(fiber_radius^2))/((x^2)) 
    count=0; 

    while abs(target_vol_fract - total_vf) > vf_tol
        if (target_vol_fract - total_vf) > 0 % Need to add fibers.
            clear fibers_old nodes_old
            message = sprintf('Creating network %i of %i. . . Adding fibers', j, num_nets);
            prog_box.Message = message;
            points_seed = points_seed + 3; % Add 5 seed points to create more fibers
            
            if strcmp(net_type, 'Delaunay') 
                [nodes_old, nodes, fibers_old, fibers] = make_aligned_dels_LB();
            elseif strcmp(net_type, 'Voronoi')
                [nodes_old, nodes, fibers_old, fibers] = make_aligned_vors(); 
            end
            
        elseif (target_vol_fract - total_vf) < 0 % Need to remove fibers.
            clear fibers_old nodes_old
            clear fibers_old nodes_old
            message = sprintf('Creating network %i of %i. . . Removing fibers', j, num_nets);
            points_seed = points_seed - 3; % Remove 3 seed points to get rid of fibers
            
            if strcmp(net_type, 'Delaunay') 
                [nodes_old, nodes, fibers_old, fibers] = make_aligned_dels_LB();
            elseif strcmp(net_type, 'Voronoi')
                [nodes_old, nodes, fibers_old, fibers] = make_aligned_vors(); 
            end
        end

        % Determine the volume fraction of the network and compare to target
        nodes1D = conv_2D_2_lin(nodes); % Reshape into 1D array
        fibers1D = conv_2D_2_lin(fibers); % Reshape into 1D array
        [lens] = calc_lens(nodes1D, fibers1D); % Returns lengths of every fiber

        int_nodes = find_int_nodes(nodes1D, boundaries); 
        num_int_nodes = length(int_nodes);
        connectivity = zeros(length(int_nodes),1);
        for j = 1:length(int_nodes)
            i = int_nodes(j);
            connectivity(i) = nnz(fibers(:) == i);
        end
        connectivity(~any(connectivity,2),:) = [];
        conn = mean(connectivity)

        total_len = sum(lens);
        total_vf = (total_len*pi*(fiber_radius^2))/((x^2)) % Box volume is 2x2x2, convert is 2.5984e8 

        count=count+1;
    end
        disp(count)
    % Plot to confirm
    plot_net(nodes, fibers);

    % Double check network orientations to see that rotation is working
    %fiber_orient_dist(nodes, fibers); % Creates orientation tensor
        [R] = calc_orient(nodes, fibers);

    length_vect(j) = total_len;


    % Write network to file
    if 1 == 0 % Text file
        fprintf('Writing info to new file... \n')
        total_nodes = size(nodes,1);
        tot_deg_freedom = 3*total_nodes;
        total_fibers = size(fibers,1);
        fnm = sprintf('%s%.2f%s%i%s', 'Del_', target_vol_fract,'_', N, '.txt');
        filename = fnm;
        fileid = fopen(filename , 'w'); % writes over any existing file

        fprintf(fileid, '%i %i %i %i %i\n', total_nodes, tot_deg_freedom, ...
            total_fibers, total_fibers, total_fibers);

        for n = 1 : total_fibers
            fprintf(fileid,'%i %i %i %f %f %f %f %f %f %f\n', ...  % file line filters
            n, ...                             % fiber num
            fibers(n,1), ...                 % node 1 num
            fibers(n,2), ...                 % node 2 num
            nodes((fibers(n,1)), 1), ...    % node 1 x coord
            nodes((fibers(n,1)), 2), ...    % node 1 y coord
            nodes((fibers(n,1)), 3), ...    % node 1 z coord
            nodes((fibers(n,2)), 1), ...    % node 2 x coord
            nodes((fibers(n,2)), 2), ...    % node 2 y coord
            nodes((fibers(n,2)), 3), ...    % node 2 z coord
            lens(n));                  % resting length of fibers (longer than distance between nodes)
        end    

        fclose(fileid);
    elseif 1 == 1 % Mat file
        fnm = sprintf('%s%.2f%s%i%s', 'Del_', target_vol_fract,'_', N, '.mat');
        save (fnm, 'nodes', 'fibers');
    end
end



%% Write text file with convert and radius values for VF
% [R] = calc_orient(nodes, fibers)
% fnm3 = sprintf('%s%.2f%s%s%s', 'AlignedDel_', target_vol_fract,'_', 'INFO', '.txt');
% filename3 = fnm3;
% fileid3 = fopen(filename3 , 'w'); % writes over any existing file
% %fprintf(fileid3,'%s %d %s %d \n', 'Convert value:', convert2, 'Minor radius:', avg_rad);
% %fprintf(fileid3,'%s %d \n', 'Convert value:', convert2);
% fclose(fileid3);

close(prog_box);

end
