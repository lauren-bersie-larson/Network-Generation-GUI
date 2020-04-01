%code to rotate fibers using jacobian matrix from <rotation.m>
%JCQ created Aug 27 2013

%n comes from n in <make_aligned_dels.m>
%nPCel=find(PC==n); %because, for example element 2206 might be the 694th
%in PC %LG edit 11/10/14 Commented this out. Don't think I need it.
for nNodes=1:size(nodes,1)
    %new_nodes(nNodes,:)=(jacobian{nPCel}*nodes(nNodes,:)')'; %use same jacobian for every node in RVE (1 jacobian/RVE)
    %LG edit 11/10/14 I'm using n instead of nPCel
    new_nodes(nNodes,:)=(jacobian{n}*nodes(nNodes,:)')';
end
clear nodes; 
%LG edit 11/10/14 comment out the next 4 lines and uncomment the 5th line
% net_axes_nodes=new_nodes;
% nodes(:,1)=net_axes_nodes(:,1);
% nodes(:,2)=(-1)*net_axes_nodes(:,3);
% nodes(:,3)=net_axes_nodes(:,2);
nodes=new_nodes;



for nNodes=1:size(nodes2,1)
    %new_nodes(nNodes,:)=(jacobian{nPCel}*nodes(nNodes,:)')'; %use same jacobian for every node in RVE (1 jacobian/RVE)
    %LG edit 11/10/14 I'm using n instead of nPCel
    new_nodes2(nNodes,:)=(jacobian{n}*nodes2(nNodes,:)')';
end
clear nodes2; 
%LG edit 11/10/14 comment out the next 4 lines and uncomment the 5th line
% net_axes_nodes=new_nodes;
% nodes(:,1)=net_axes_nodes(:,1);
% nodes(:,2)=(-1)*net_axes_nodes(:,3);
% nodes(:,3)=net_axes_nodes(:,2);
nodes2=new_nodes2;
