% Removes free fibers not connected to any other fiber within the network

[n_nodes, ~] = size(nodes);
[n_fibers, ~] = size(fibers);



node_count = zeros(n_nodes,1);

for i = 1:n_fibers
   
    node1 = fibers(i,1);
    node2 = fibers(i,2);
    node_count(node1,1) = node_count(node1,1) + 1;
    node_count(node2,1) = node_count(node2,1) + 1;
    
end

for i = 1:n_fibers
   
    node1 = fibers(i,1);
    node2 = fibers(i,2);
    
    
    if (node_count(node1,1) == 1 && node_count(node2,1) == 1)
        
        nodes(node1,3) = 0;
        node_count(node1,1) = 2;
        
    end
        
        
    
end


% new_fibers = zeros(n_fibers-count, 2);
% j = 1;
% for i = 1:n_fibers
%     
%    
%     if fib_reg(i,1) == 1
%         
%         new_fibers(j,:) = fibers(i,:);
%         j = j + 1;
%     end
%     
% end
