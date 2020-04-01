function plot_net1(nodes, fibers, fib_reg)

%file = 'C:\Users\Lauren\Google Drive\Lab Work\Kidney Project\New Kidney Code\Const Flow\Flow_4e-17\\Healthy 0.08 Networks\DELS\Aligned_net_1.txt';
%[nodes, fibers] = read_network(file);
% plot_net(nodes, fibers)
%
% in:
%
% nodes            N x 3 coordinates for N nodes
% fibers           N x 2 start-end nodes for N fibers
%
% last rev:
%
% Mon 2 June 2014 RD


%figure;

plot3(nodes(:,1) , nodes(:,2), nodes(:,3),'o', 'LineWidth',5, ...
    'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1);

%axis_val = max( nodes(1:end) ) + 0.1;     % find a max/min for plot
%axis([-axis_val axis_val -axis_val axis_val -axis_val axis_val]);
%axis([-1 +1 -1 +1 -1 +1]);

%hold on;

for n = 1 : size(fibers, 1) % count rows
    
    x(1) = nodes(fibers(n,1), 1); % node 1 x coord
    y(1) = nodes(fibers(n,1), 2); % node 1 y coord
    z(1) = nodes(fibers(n,1), 3); % node 1 z coord
    
    x(2) = nodes(fibers(n,2), 1); % node 2 x coord
    y(2) = nodes(fibers(n,2), 2); % node 2 y coord
    z(2) = nodes(fibers(n,2), 3); % node 2 z coord
    
% %     if fib_reg(n) == 1
% %         plot3(x, y, z, 'k','LineWidth',2);
% %     elseif fib_reg(n) == 2
         %color = [185 211 238]/255;
         plot3(x, y, z, 'k','LineWidth', .5);
% %     elseif fib_reg(n) == 3
% %         plot3(x, y, z, 'g', 'LineWidth', 1);
% %     elseif fib_reg(n) == 0
% %         % Do nothing
% %     end
    
    hold on;
    
end

set(gcf, 'color', 'white');
axis equal;
set(gca,'xtick',[],'ytick',[], 'ztick', [])
%axis( [-1 +1 -1 +1 -1 +1] );

%axis_lim = max( nodes(1:end) );

%axis( [-axis_lim +axis_lim -axis_lim +axis_lim -axis_lim +axis_lim] );

%axis([0 2 0 2 0 2 ]);



%xlabel('x'); ylabel('y'); zlabel('z');

end
