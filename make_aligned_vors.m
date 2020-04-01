function [nodes_old, nodes, fibers_old, fibers] = make_aligned_vors()
% make many aligned dels

%CHECK
%that saving networks to the right directory if you want to save them to a directory
%if saving images, check that they're not erasing other important images

% Updated 7-20-17 to include global points_seed parameter to adjust volume
% fraction, LMB

global points_seed 
%%
filename = fullfile('net_parameters.txt');
fid_parameters = fopen(filename, 'r');
parameters = fscanf(fid_parameters, ['%s %f']);

lambdax = parameters(9);
lambday = parameters(18);
lambdaz = parameters(27);

xmin = parameters(32);
xmax = parameters(37);
ymin = parameters(42);
ymax = parameters(47);
zmin = parameters(52);
zmax = parameters(57);
%N = parameters(59);

N = 1;

for n = 1:N %to make N networks change to: n=1:N
    pts_xyz =(2.0 .* rand(points_seed, 3) - 1.0); % N random pts from -1 to +1
    % pts_xyz = 2.0 .* rand(2000, 3) - 1.0; % N random pts from -1 to +1

    [nodes_old, fibers_old] = make_vor(pts_xyz); % fnxn in netmat


    % geometrically stretch fibers in x/y/z
    nodes_old(:,1)=nodes_old(:,1)*lambdax; %nodes is N x 3 for N nodes
    nodes_old(:,2)=nodes_old(:,2)*lambday;
    nodes_old(:,3)=nodes_old(:,3)*lambdaz;
    
    clipbox = [xmin, xmax, ymin, ymax, zmin, zmax]; % xmin xmax ymin ymax zmin zmax
    [nodes, fibers] = clip_net(nodes_old, fibers_old, clipbox); % fnxn in netmat
  
    % get giant network
    [nodes, fibers] = get_giant(nodes, fibers);

    fname = sprintf('net_%i.del', n);
    %fname = sprintf('network_library/net_%i.del', n); %use if saving many nets to network_library directory

    fprintf('creating %s\n\n', fname);
    %put_net(nodes, fibers, fname, ones(length(fibers),1));% This writes out the network file
  
    % calculate network fiber orientation

    [R] = calc_orient(nodes, fibers);

    disp( size(fibers,1) ) %number of fibers in network
      
    plot_net1(nodes,fibers, ones(length(fibers))); %display network

    %%SAVE NETWORK(uncomment this section below if want to save .jpg and/or .fig  images of network)   
    %     fpath='images/';
  
    %     filename_jpg=sprintf('net_%d.jpg',n); %save jpg image
    %     saveas(gca,fullfile(fpath,filename_jpg));

    %     filename_fig=sprintf('net_%d.fig',n); %save fig image
    %     saveas(gca,fullfile(fpath,filename_fig));

end

%rmpath('netmat/');
end