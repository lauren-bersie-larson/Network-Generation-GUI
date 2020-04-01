%Aug 15, 2013

%Extract center of element

%CHECK:
%The data directory is correct folder
clc; clear all; close all;

%DATA_DIR = '../../Matlab plots/';  %check if this should be PC_GEO or PC_GEO_OPEN

nop_file = ['nop'];
nodes_file = ['nodes'];
nop = dlmread(nop_file, ' ', 1, 0);
nop = nop(:, 2:9); 
nodes= dlmread(nodes_file, ' ', 1, 1); 
num_els = size(nop, 1);
old_nodes=nodes;

%Change coordinates from code [1,2,3] to Prolate Spheroid [x,y,z]
%z must change signs (see page 20 of lab notebook)
%center(:,2) is their z, which neds to become negative
x_coor=nodes(:,1);
y_coor=nodes(:,2); %LG edit 11/11/14 JQ had 3 instead of 2 fo r the column number
%z_coor=(-1).*nodes(:,2); %LG edit 11/10/14
z_coor = nodes(:,3);
clear nodes;
nodes=horzcat(x_coor,y_coor,z_coor);


%a=.4; %prolate spheroidal constant
%a=sqrt(.5^2-.25^2); %LG edit 11/10/14

jacobian=cell(num_els,1,1);
%Find center of element
for nEl=1:num_els
    nodes_el=nop(nEl,:);
    loc_nodes=nodes(nodes_el',:);
    center(nEl,:)=mean(loc_nodes,1);
    if (nEl < 8)
    center(nEl,1:2)=center(nEl,1:2)-10; %LG edit March 20, 2015. To make the (0,0) point at the center of the capillary
    elseif ((nEl>7) && (nEl<15))
        center(nEl,1:2)=center(nEl,1:2)-10.3; %LG edit March 20, 2015. To make the (0,0) point at the center of the capillary
    else
        center(nEl,1:2)=center(nEl,1:2)-10.6; %LG edit March 20, 2015. To make the (0,0) point at the center of the capillary
    end
    lengthh=1; %<<<<<<<<<<-------------- don't understand why this is 1 or what it's doing, but it might be unused LG edit 11/10/14
    %Make unit vector
    %length(nEl)=sqrt(center(nEl,1)^2+center(nEl,2)^2+center(nEl,3)^2);
    %length(nEl)=sqrt(center(nEl,1)^2+center(nEl,2)^2); %LG edit 11/10/14 %Julia said that I can probably ignore it. If something doesn't work, try setting it to 1.
end
%%
%For practice
% %length=(sqrt(1^2+.5^2+.5^2));
 %length=1; %LG edit 11/10/14
% % %center=[-.25 0 0;.25 0 0; 0 0 0; 0 0 .25; 0 .5 0; 0 -.5 0; -.25 0 0; -.1 0 0];
% % %center=[-.25 .0001 0.0001; .25 0.0001 0.0001;.0001 0.0001 0.0001; 0.0001 .25 0.0001; 0.0001 0.0001 -.5; 0.0001 0.0001 .5; -.25 0.0001 0.0001; -.1 0.0001 0.0001; .125 .0001 .4 ;.125 .0001 -.4];
% %center=[.25 0 0; .1768 0 0.3536; .1768 0 -.3536; .1768 .1768 0; -.1768 .1768 0; 0 .25 0];
% center=[0 .25 0; 0 .176775 .353552; 0 .176775 -.353552; 0 .0005 0.4999];
% num_els=size(center,1);
%%
for nEl=1:num_els
    center_u(nEl,:)=center(nEl,:)./max(lengthh); %is Julia doing this so all numbers will be between 0 and 1? So then I should use the length calcualted in line 39? LG comment 11/10/14
    %Calculate Jacobian
    %sigma(nEl)=(1/(2*a))*(sqrt(center_u(nEl,1)^2+center_u(nEl,2)^2+(center_u(nEl,3)+a)^2)+sqrt(center_u(nEl,1)^2+center_u(nEl,2)^2+(center_u(nEl,3)-a)^2));
    %LG edit 11/10/14
    %tau(nEl)=(1/(2*a))*(sqrt(center_u(nEl,1)^2+center_u(nEl,2)^2+(center_u(nEl,3)+a)^2)-sqrt(center_u(nEl,1)^2+center_u(nEl,2)^2+(center_u(nEl,3)-a)^2));
    %LG edit 11/10/14
    %phi(nEl)=atan2(center_u(nEl,2),center_u(nEl,1)); %used atan2 instead of atan %LG edit 11/10/14
    theta(nEl) = atan(center_u(nEl,2)/center_u(nEl,1)); %LG edit 11/10/14
    th_deg=theta(nEl)*180/pi;
    %mu(nEl)=acosh(sigma(nEl)); %LG edit 11/10/14
    %nu(nEl)=acos(tau(nEl)); %LG edit 11/10/14
    r(nEl) = sqrt(center_u(nEl,1)^2+center_u(nEl,2)^2); %LG_edit 11/10/14 should be the same everywhere. Just the radius. Normalized?
%     if center_u(nEl,2)==0 & center_u(nEl,1)>=0
%         phi(nEl)=0;
%     end
    if center_u(nEl,1)==0
        theta(nEl)=pi/2; %LG edit 11/10/14. This used to be phi, changed it to my theta
    end
    %Create jacobian values
    %LG's equations
    xu(nEl) = cos(theta(nEl));
    yu(nEl) = sin(theta(nEl));
    zu(nEl) =0;
    xv(nEl) = -r(nEl)*(sin(theta(nEl)));
    yv(nEl) = r(nEl)*cos(theta(nEl));
    zv(nEl) = 0;
    xz(nEl) = 0;
    yz(nEl) = 0;
    zz(nEl) = 1;
    
    normu(nEl) = sqrt(xu(nEl)^2+yu(nEl)^2+zu(nEl)^2);
    normv(nEl) = sqrt(xv(nEl)^2+yv(nEl)^2+zv(nEl)^2);
    normz(nEl) = sqrt(xz(nEl)^2+yz(nEl)^2+zz(nEl)^2);
    
    %comment out JQ's equaitons
%         xu(nEl)=a*cosh(mu(nEl))*sin(nu(nEl))*cos(phi(nEl));
%         yu(nEl)=a*cosh(mu(nEl))*sin(nu(nEl))*sin(phi(nEl));
%         zu(nEl)=a*sinh(mu(nEl))*cos(nu(nEl));
%         xv(nEl)=a*sinh(mu(nEl))*cos(nu(nEl))*cos(phi(nEl));
%         yv(nEl)=a*sinh(mu(nEl))*cos(nu(nEl))*sin(phi(nEl));
%         zv(nEl)=-a*cosh(mu(nEl))*sin(nu(nEl));
%         xp(nEl)=-a*sinh(mu(nEl))*sin(nu(nEl))*sin(phi(nEl));
%         yp(nEl)=a*sinh(mu(nEl))*sin(nu(nEl))*cos(phi(nEl));
%         zp(nEl)=0;
%         normu(nEl)=sqrt(xu(nEl)^2+yu(nEl)^2+zu(nEl)^2);
%         normv(nEl)=sqrt(xv(nEl)^2+yv(nEl)^2+zv(nEl)^2);
%         normp(nEl)=sqrt(xp(nEl)^2+yp(nEl)^2+zp(nEl)^2);
        %jacobian{nEl}=[[xu(nEl) xv(nEl) xp(nEl)]./norm(nEl,1); [yu(nEl) yv(nEl) yp(nEl)]./norm(nEl,2); [zu(nEl) zv(nEl) zp(nEl)]./norm(nEl,3)];
        %jacobian{nEl}=[xu(nEl)/normu(nEl) xv(nEl)/normv(nEl) xp(nEl)/normp(nEl);yu(nEl)/normu(nEl) yv(nEl)/normv(nEl) yp(nEl)/normp(nEl); zu(nEl)/normu(nEl) zv(nEl)/normv(nEl) zp(nEl)/normp(nEl)];
        %LG edit 11/10/14 JQ had the above eqn. Mine is 2x2 below
        jacobian{nEl} = [xu(nEl)/normu(nEl) xv(nEl)/normv(nEl) xz(nEl)/normz(nEl); yu(nEl)/normu(nEl) yv(nEl)/normv(nEl) yz(nEl)/normz(nEl); zu(nEl)/normu(nEl) zv(nEl)/normv(nEl) zz(nEl)/normz(nEl)];
        
        %jacobian{nEl}=[xu(nEl) xv(nEl) xp(nEl); yu(nEl) yv(nEl) yp(nEl); zu(nEl) zv(nEl) zp(nEl)];
end


    % %Dot product
% for nEl=1:num_els
%     cos_alpha(nEl)=dot(center_u(nEl,:)',[1;0;0]);
%     cos_beta(nEl)=dot(center_u(nEl,:)',[0;1;0]);
%     cos_gamma(nEl)=dot(center_u(nEl,:)',[0;0;1]);
% end