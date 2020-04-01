% Removes free fibers not connected to any other fiber within the network

[n_nodes, ~] = size(nodes);
[n_fibers, ~] = size(new_fibers);



node_count = zeros(n_nodes,1);

for i = 1:n_fibers
   
    node1 = new_fibers(i,1);
    node2 = new_fibers(i,2);
    node_count(node1,1) = node_count(node1,1) + 1;
    node_count(node2,1) = node_count(node2,1) + 1;
    
end




% for i = 1:n_fibers
%     
%     node1 = new_fibers(i,1);
%     node2 = new_fibers(i,2);
%     
%     if (node_count(node1,1) == 1 && node_count(node2,1) == 1)
%         
%         node_count(node1,1) = 0;
%         node_count(node2,1) = 0;
%         
%     end
%     
% end

j = 1;
for i = 1:n_nodes
    
    
    if ~(node_count(i,1) == 0)
        
        new_nodes(j,:) = nodes(i,:);
        node_numbers(j,1) = i;
        j = j + 1;
        
    end
        
        
    
end

[new_n,~] = size(node_numbers);

for i = 1:new_n
    
    indices = find(new_fibers == node_numbers(i));
    new_fibers(indices) = i;
    
end
