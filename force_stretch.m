% STRETCH A MICROSCALE NETWORK
%
% -deform the boundaries
% -solve for the nodal displacements <vect> that gives us F=0 on the nodes
%
% FIXED POISSON RATIO OF ZERO -- FIXED FACE STRETCH


function [nodes_n,fibers_n,S11,S22,S33,stretch, face_forces]=force_stretch(nodes,fibers,fib_type)

for n = 1:1 % stretch a single network
    
    [nodes_final,fibers_final, S11,S22,S33,stretch, fib_reg, face_forces]=stretch_net(n,nodes,fibers, fib_type);
    
end

% Section to plot final network
if 1 == 1
    % Convert nodes and fibers to be matrices again, not 1D vectors
    
    nodes_n = zeros((length(nodes_final)/3),3);
    fibers_n = zeros((length(fibers_final)/2),2);
    
    nodes_n(:,1) = nodes_final(1:3:end);
    nodes_n(:,2) = nodes_final(2:3:end);
    nodes_n(:,3) = nodes_final(3:3:end);
    
    fibers_n(:,1) = fibers_final(1:2:end);
    fibers_n(:,2) = fibers_final(2:2:end);
    % Plot final network
%     figure;
%     plot_net(nodes_n, fibers_n, fib_reg);
end

S11 = [0,S11];
S22 = [0,S22];
S33 = [0,S33];
stretch = [1,stretch];


end



function [nodes,fibers,S_11,S_22,S_33,M_GS, fib_reg, face_forces] = stretch_net(iter_num,nodes,fibers, fib_type)


stretch_steps = 40;
fib_r = 100e-9;
fib_limit = [1.35,1.8,1.375, 100];

% Different cut off parameters
% 0 is off, 1 is on

strain_cutoff = 1;
stress_cutoff = 0;

% stretch net
% Constant displacement steps
s = 1.01.*ones(1,stretch_steps); % 2 percent stretch over N steps
% s(1:40) = 1.02;

% Vary displacement sinusoidally
% angle = linspace(0, 4*pi, stretch_steps);
%
% s = 1.0 + 0.1*sin(angle);

  
% convert 2D ==> 1D arrays

nodes = conv_2D_2_lin(nodes);

fibers = conv_2D_2_lin(fibers);



boundaries = [-0.5 +0.5 -0.5 +0.5 -0.5 +0.5];
% boundaries = [-1 +1 -1 +1 -1 1]; % x y z

bnd_node_nums = get_bnd_nodes(nodes, boundaries);



init_lens = calc_lens(nodes, fibers);



num_fibers = length(fibers) / 2;



fib_mods = ones(1, num_fibers);    % modulus for each fiber
fib_bs   = ones(1, num_fibers);    % B for each fiber



% fiber modulus and B values
for i=1:num_fibers
    if fib_type(i) == 1
        fib_mods(i) = 0.15e5;
        fib_bs(i) = 0.1;
    elseif fib_type(i) == 2
        fib_mods(i) = 0.35e6;
        fib_bs(i) = 1.4;
    elseif fib_type(i) == 3
        fib_mods(i) = 1.74e6;
        fib_bs(i) = 2.6;
    elseif fib_type(i) == 4
        fib_mods(i) = 0.15e5; %7.4;
        fib_bs(i) = 0.1; %0.01;
    end
end


total_length = sum(init_lens);


fib_areas =  pi .* fib_r .* fib_r .* ones(1, num_fibers); %


Fiber_volume = total_length*pi*fib_r*fib_r;

Box_volume = (boundaries(2) - boundaries(1))*(boundaries(4) - boundaries(3))*...
    (boundaries(6) - boundaries(5));

fiber_vol_fract = 0.05; %Fiber_volume/Box_volume;

prot_conc = 2.0; % mg/ml

count = 0;
count_c = 0;
count_e = 0;
count_ilc = 0;

fib_reg = ones(1,length(fibers));

face_forces = zeros(length(s),1);

for n = 1:length(s)
    

    rve_stretch(1) = 1/sqrt(s(n));              % in x -- slightly different
    rve_stretch(2) = s(n);                 % in y -- no poisson   
    rve_stretch(3) = 1/sqrt(s(n));                 % in z -- no poisson
    
    [nodes, fibers, S1, face_forces(n), Stress, fib_stress] = solve_fixed_cauchy(nodes, fibers, bnd_node_nums, ...
        init_lens, rve_stretch, fib_mods, fib_areas, fib_bs, ...
        prot_conc, total_length, fiber_vol_fract, 2);
    

    % CAUCHY STRESS
    
    S_11(n) = S1(1,1);
    S_22(n) = S1(2,2);
    S_33(n) = S1(3,3); 
    
    
    % Fiber breaking part
    % Strain breakage
    
    
    if strain_cutoff == 1
        new_length = calc_lens(nodes, fibers);
        fib_stretch = new_length./init_lens;
        
        for m = 1:length(fib_stretch)
            
            if fib_stretch(m) > fib_limit(fib_type(m))
                
                fib_mods(m) = 2.51e-6;
                fib_reg(m) = 0;
                count = count + 1;
                
                switch (fib_type(m))
                    case 1
                        count_c = count_c + 1;
                    case 2
                        count_e = count_e + 1;
                    case 3
                        count_ilc = count_ilc + 1;
                end
                
            end
        end
    end
    if stress_cutoff == 1
        
    end
    
    if (strain_cutoff == 1 || stress_cutoff == 1)
        
        fprintf('Run : %d  Fibers broken : %d\n', n,count);
        fprintf('Max lambda : %f\n\n\n', max(fib_stretch));
    
    else
        fprintf('Run : %d\n\n', n);
        fprintf('Max fib stress: %3.2e\n\n', max(fib_stress));
    end
    
    % Check to see if required number of fibers are broken
      if count > 250
          break;
      end
    
    % Check to see if stress drops suddenly
%     if ((n > 10) && (count > 1))
%         if (S_22(n) > S_22(n-3))
%             fprintf('Stopped due to stress drop\n');
%             fprintf('Fibers broken:\n');
%             fprintf('C : %f, E : %f, ILC : %f\n', 100*count_c/count...
%                 , 100*count_e/count, 100*count_ilc/count);
%             break;
%         end
%     end
    
end

% S_11 = S_11 - S_33;%subtract z-stresses from all (hydrostatic pressure term)
% S_22 = S_22 - S_33;

M_GS = cumprod(s); % MODEL LAMBDA VECTOR


% figure;
% hold on;
% 
% plot(M_GS(1:length(S_11)), S_11, '-xk');
% plot(M_GS(1:length(S_22)), S_22, '-xm');
% plot(M_GS(1:length(S_33)), S_33, '-xg');





end
