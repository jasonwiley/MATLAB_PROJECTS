clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT (READ THIS SECTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Simulation box ---
DSx = 14.4;
DSy = 15.8;
DSz = 13.20857;


Lx = DSx*3;
Ly = DSy*3;
Lz = DSz*3;
box = [Lx, Ly, Lz];

% --- Lamella parameters ---
nLamellae = 3;             % number of lamellae
domain_spacing = DSz;      % spacing between lamella centers
lamella_axis = 'z';         % orientation axis
a_peak = 7;            % center of first peak
a_min  = 6.135;         % min spacing
a_max  = 8.385; 
% Define number of particles per lamella individually
N_per_lamella = [39, 39, 39];  % must sum to total N
positions = [];
dx = 0;  % example shift in X
dy = 4; % example shift in Y
dz = 0;  % example shift in Z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Compute lamella centers ---
start = -domain_spacing*(nLamellae-1)/2;
lamella_centers = start + (0:nLamellae-1)*domain_spacing;

% --- Generate particles per lamella ---
for k = 1:nLamellae
    center = lamella_centers(k);
    N_this_lamella = N_per_lamella(k);
    
    % Determine plane axes based on lamella orientation
    switch lamella_axis
        case 'x'
            plane_L1 = Ly; plane_L2 = Lz;  % plane: Y-Z
        case 'y'
            plane_L1 = Lx; plane_L2 = Lz;  % plane: X-Z
        case 'z'
            plane_L1 = Lx; plane_L2 = Ly;  % plane: X-Y
        otherwise
            error('Invalid lamella axis. Choose x, y, or z.');
    end
    
    % Estimate number of points along plane axes
    n1 = ceil(sqrt(N_this_lamella*plane_L1/plane_L2));
    n2 = ceil(N_this_lamella/n1);
    
    d1 = plane_L1 / n1;
    d2 = plane_L2 / n2 * sqrt(3)/2;
    
    [G1, G2] = meshgrid(0:n1-1, 0:n2-1);
    G1(2:2:end,:) = G1(2:2:end,:) + 0.5;
    G1 = G1 * d1 - plane_L1/2 + d1/2;
    G2 = G2 * d2 - plane_L2/2 + d2/2;
    
    hex_pos = [G1(:), G2(:)];
    hex_pos = hex_pos(1:min(N_this_lamella, size(hex_pos,1)),:);
    % --- Randomize in-plane spacing to match g(r) peak 6.135-8.385 Å ---
            % max spacing
    
    n_particles = size(hex_pos,1);
    scale_factors = a_min/a_peak + (a_max/a_peak - a_min/a_peak)*rand(n_particles,1);
    hex_pos = hex_pos .* scale_factors;
    % Map plane coordinates to 3D
    switch lamella_axis
        case 'x'
            Y = hex_pos(:,1);
            Z = hex_pos(:,2);
            X = center*ones(size(Y));
        case 'y'
            X = hex_pos(:,1);
            Z = hex_pos(:,2);
            Y = center*ones(size(X));
        case 'z'
            X = hex_pos(:,1);
            Y = hex_pos(:,2);
            Z = center*ones(size(X));
    end
    
    lamella_pos = [X, Y, Z];
    positions = [positions; lamella_pos];
end

fprintf('Generated %d particles in %d lamellae along %s-axis.\n', size(positions,1), nLamellae, lamella_axis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT (READ THIS SECTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Add Gaussian noise ---
sigmaX = 0.23; 
sigmaY = 0.5;
sigmaZ = 0.23;
positions(:,1) = positions(:,1) + sigmaX * randn(size(positions,1),1);
positions(:,2) = positions(:,2) + sigmaY * randn(size(positions,1),1);
positions(:,3) = positions(:,3) + sigmaZ * randn(size(positions,1),1);

%--- Introduce defects ---
% defect1 = 8; 
defect2 = 55; 
% defect3 = 62;
% positions(defect1,3) = positions(defect1,3) + 3;
positions(defect2,2) = positions(defect2,2) + 3;
% positions(defect3,3) = positions(defect3,3) -7;
% Ensure all particles stay inside the box


positions = positions + [dx, dy, dz];
positions(:,1) = min(max(positions(:,1), -Lx/2 + 0.2), Lx/2 - 0.2);
positions(:,2) = min(max(positions(:,2), -Ly/2 + 0.2), Ly/2 - 0.2);
positions(:,3) = min(max(positions(:,3), -Lz/2 + 0.2), Lz/2 - 0.2);


% --- Save ---
save('F_hex_lamellae_positions.mat','positions');
writematrix(positions,'F_hex_lamellae_positions.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_boundary = [-Lx/2 Lx/2];
Y_boundary = [-Ly/2 Ly/2];
Z_boundary = [-Lz/2 Lz/2];
% --- Plot particles and box ---
figure;
scatter3(positions(:,1), positions(:,2), positions(:,3), 60, 'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('%d Lamellae along %s-axis (hexagonal lattice) with box', nLamellae, lamella_axis));
axis equal; grid on;
hold on;

% Box corners
corners = [
    -Lx/2, -Ly/2, -Lz/2;
     Lx/2, -Ly/2, -Lz/2;
     Lx/2,  Ly/2, -Lz/2;
    -Lx/2,  Ly/2, -Lz/2;
    -Lx/2, -Ly/2,  Lz/2;
     Lx/2, -Ly/2,  Lz/2;
     Lx/2,  Ly/2,  Lz/2;
    -Lx/2,  Ly/2,  Lz/2
];

% Box edges
edges = [
    1 2; 2 3; 3 4; 4 1;   % bottom face
    5 6; 6 7; 7 8; 8 5;   % top face
    1 5; 2 6; 3 7; 4 8    % vertical edges
];

for e = 1:size(edges,1)
    plot3(corners(edges(e,:),1), corners(edges(e,:),2), corners(edges(e,:),3), 'k-', 'LineWidth', 1.5);
end

hold off;
