[rows, ~] = size(nodes_old);

for i = 1:rows
    
    if (nodes_old(i,3) < 0.7 && nodes_old(i,3) > -0.7)
        
        nodes_old(i,3) = 0;
        
    else
        if nodes_old(i,3) > 0
            nodes_old(i,3) = nodes_old(i,3) - 1;
        else
            nodes_old(i,3) = nodes_old(i,3) + 1;
        end
        
    end
    
end