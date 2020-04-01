function [stress] = calc_net_str(nodes, forces, bnd_node_nums)

% bnd_node_nums -- 1xN array of node numbers for N boundary nodes
% forces -- 1x3N array of xyzxyz... force components acting on N nodes
% nodes -- 1x3N array of xyzxyz... coordinates for N nodes
%
% stress -- 3x3 array for the volume averaged stress tensor

% Remove nodes with really high forces from calcualtion
% ind1=find(bnd_node_nums==121);
% ind2=find(bnd_node_nums==125);
% bnd_node_nums(ind2)=[];
% bnd_node_nums(ind1)=[];

% Sum xi * fi through matrix multiplication
xfx = nodes(3*bnd_node_nums-2) * forces(3*bnd_node_nums-2)';
yfx = nodes(3*bnd_node_nums-1) * forces(3*bnd_node_nums-2)';
zfx = nodes(3*bnd_node_nums-0) * forces(3*bnd_node_nums-2)';

xfy = nodes(3*bnd_node_nums-2) * forces(3*bnd_node_nums-1)';
yfy = nodes(3*bnd_node_nums-1) * forces(3*bnd_node_nums-1)';
zfy = nodes(3*bnd_node_nums-0) * forces(3*bnd_node_nums-1)';

xfz = nodes(3*bnd_node_nums-2) * forces(3*bnd_node_nums-0)';
yfz = nodes(3*bnd_node_nums-1) * forces(3*bnd_node_nums-0)';
zfz = nodes(3*bnd_node_nums-0) * forces(3*bnd_node_nums-0)';

%Make it incomperssible
% xfx=xfx-zfz;
% yfy=yfy-zfz;
% zfz=zfz-zfz;

%Incompressible for uniaxial
%xfx=xfx-yfy;
%check if the y and z forces are the same. They should be


% Formula for this is Sij=1/V*Sum(xiFj); The off diagnonal terms are
% symmetric, e.g. xfy=yfx. However, there is a very small difference in
% these numbers due to numerical issues, so we just average them. The
% difference seems to be at the 10-14 decimal place, so probably could
% ignore.

% compact notation [S11 S12 S13 S22 S23 S33]
%stress= [xfx; 0.5*(xfy+yfx); 0.5*(zfx+xfz); yfy;0.5*(yfz+zfy); zfz]; 

% full notation [S11 S12 S13; S21 S22 S23; S31 S32 S33]
stress= [ xfx 0.5*(xfy+yfx) 0.5*(zfx+xfz); 
          0.5*(xfy+yfx) yfy 0.5*(yfz+zfy); 
          0.5*(zfx+xfz) 0.5*(yfz+zfy) zfz ];

end
