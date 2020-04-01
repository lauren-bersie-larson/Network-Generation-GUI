function plot_net(nodes, fibers)


% plot_net(nodes, fibers)
%
% in:
%
% nodes            N x 3 coordinates for N nodes
% fibers           N x 2 start-end nodes for N fibers
%
% last rev:
%
% tue nov 6 2012 mfh


figure;

plot3(nodes(:,1) , nodes(:,2), nodes(:,3),'o', 'LineWidth',0.2, ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.8],'MarkerSize',3);

%axis_val = max( nodes(1:end) ) + 0.1;     % find a max/min for plot
%axis([-axis_val axis_val -axis_val axis_val -axis_val axis_val]);
%axis([-1 +1 -1 +1 -1 +1]);

hold on;

for n = 1 : size(fibers, 1) % count rows
    
    x(1) = nodes(fibers(n,1), 1); % node 1 x coord
    y(1) = nodes(fibers(n,1), 2); % node 1 y coord
    z(1) = nodes(fibers(n,1), 3); % node 1 z coord
    
    x(2) = nodes(fibers(n,2), 1); % node 2 x coord
    y(2) = nodes(fibers(n,2), 2); % node 2 y coord
    z(2) = nodes(fibers(n,2), 3); % node 2 z coord
    
    plot3(x, y, z, 'b');
    
    hold on;
    
end

set(gcf, 'color', 'white');
axis equal;
%axis( [-1 +1 -1 +1 -1 +1] );

axis_lim = max( nodes(1:end) );

axis( [-axis_lim +axis_lim -axis_lim +axis_lim -axis_lim +axis_lim] );

%axis([0 2 0 2 0 2 ]);



xlabel('n1')'; ylabel('n2'); zlabel('n3');

end
