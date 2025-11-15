clear all
close all

system = load('system.mat');

A_index = find(system.attype == 'A');
B_index = find(system.attype == 'B');
C_index = find(system.attype == 'C');
D_index = find(system.attype == 'D');

A_pos = system.pos(A_index, :);
B_pos = system.pos(B_index, :);
C_pos = system.pos(C_index, :);
D_pos = system.pos(D_index, :);
% A beads positions (possibly wrapped)
% A_index = find(system.attype == 'A');
% A_pos = system.pos(A_index,:);       % N x 3
% box = system.dim;                    % 1 x 3
% 
% % Initialize unwrapped positions
% A_unwrapped = zeros(size(A_pos));
% 
% % Anchor first bead
% A_unwrapped(1,:) = A_pos(1,:);
% 
% % Loop over remaining beads and cumulatively unwrap
% for i = 2:size(A_pos,1)
%     % Compute 3D displacement vector to previous bead
%     dr = A_pos(i,:) - A_pos(i-1,:);
% 
%     % Apply full 3D minimum-image correction
%     dr = dr - box .* round(dr ./ box);
% 
%     % Cumulatively reconstruct unwrapped position
%     A_unwrapped(i,:) = A_unwrapped(i-1,:) + dr;
% end
% A_pos=A_unwrapped;
box = system.dim;   % [Lx Ly Lz], already 1x3
% 
% % % --- Apply PBC minimum-image to all positions ---
A_pos = A_pos - box .* round(A_pos ./ box);
B_pos = B_pos - box .* round(B_pos ./ box);
C_pos = C_pos - box .* round(C_pos ./ box);
D_pos = D_pos - box .* round(D_pos ./ box);
% 

center_vector = sqrt(sum(A_pos.^2, 2));
[~, index] = min(center_vector);
A_pos = A_pos - A_pos(index, :);

% Parameters
lower_bound = -5;
upper_bound = 5;
sensitivity = 1;
angles = 0:sensitivity:180;
num_angles = numel(angles);

% Preallocate
std_vals = nan(num_angles, num_angles);

for i = 1:num_angles     % X rotation
    rotated_x = A_pos * rotx(angles(i))';
    
    for k = 1:num_angles % Z rotation
        rotated_x_z = rotated_x * rotz(angles(k))';
        
        % Points in slice
        in_idx = rotated_x_z(:,2) >= lower_bound & rotated_x_z(:,2) <= upper_bound;
        
        if any(in_idx)
            std_vals(i,k) = std(rotated_x_z(in_idx,2), 'omitnan');
        end
    end
end

% Find the minimum std (ignore NaNs)
[min_val, linear_idx] = min(std_vals(:), [], 'omitnan');
[ix, iz] = ind2sub(size(std_vals), linear_idx);

% Map indices to angles
best_angle_x = angles(ix);
best_angle_z = angles(iz);

fprintf('Best angles: X = %.1f°, Z = %.1f°, min std = %.4f\n', ...
    best_angle_x, best_angle_z, min_val);

% Apply best rotation
rotated_A = A_pos * rotx(best_angle_x)' * rotz(best_angle_z)';
rotated_C = C_pos * rotx(best_angle_x)' * rotz(best_angle_z)';
rotated_B = B_pos * rotx(best_angle_x)' * rotz(best_angle_z)';
system.pos = system.pos * rotx(best_angle_x)' * rotz(best_angle_z)';
% figure()
% scatter3(rotated_A(:,1), rotated_A(:,2), rotated_A(:,3))
% rotated_D = D_pos * rotx(best_angle_x)' * rotz(best_angle_z)';
% Ly = system.dim(2);
% y = rotated_A(:,2);
% 
% % % shift y > Ly/2 down by Ly
% y(y > Ly/2) = y(y > Ly/2) - Ly;
% 
% % shift y < -Ly/2 up by Ly
% y(y < -Ly/2) = y(y < -Ly/2) + Ly;
% 
% % put back into rotated_A
% rotated_A(:,2) = y;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optflag=2;
first_lamellae_lower=-20;
first_lamellae_upper=-10;
second_lamellae_lower=-5;
second_lamellae_upper=5;
if optflag==0
    optional_third_lower = 5;
    optional_third_upper = 20;
elseif optflag==2
    optional_third_lower = 10;
    optional_third_upper = 20;
    optional_fourth_lower = 20;
    optional_fourth_upper = 40;
end
if optflag == 0
    first_lamellae  = rotated_A(:,2) >= first_lamellae_lower  & rotated_A(:,2) <= first_lamellae_upper;
    second_lamellae = rotated_A(:,2) >= second_lamellae_lower & rotated_A(:,2) <= second_lamellae_upper;
    first_center  = mean(rotated_A(first_lamellae, 2));
    second_center = mean(rotated_A(second_lamellae, 2));
    domain_spacing = (abs(second_center - first_center));
    plane_SD=mean([std(rotated_A(first_lamellae, 2)),std(rotated_A(second_lamellae, 2))]);
elseif optflag==1
    first_lamellae  = rotated_A(:,2) >= first_lamellae_lower  & rotated_A(:,2) <= first_lamellae_upper;
    second_lamellae = rotated_A(:,2) >= second_lamellae_lower & rotated_A(:,2) <= second_lamellae_upper;
    first_center  = mean(rotated_A(first_lamellae, 2));
    second_center = mean(rotated_A(second_lamellae, 2));
    optional_third_lamellae = rotated_A(:,2) >= optional_third_lower & rotated_A(:,2) <= optional_third_upper;
    third_center = mean(rotated_A(optional_third_lamellae, 2));
    domain_spacing = (abs(third_center - first_center) + abs(second_center-first_center))/2;
    plane_SD=mean([std(rotated_A(first_lamellae, 2)),std(rotated_A(second_lamellae, 2)), std(rotated_A(optional_third_lamellae,2))]);
elseif optflag==2
    first_lamellae  = rotated_A(:,2) >= first_lamellae_lower  & rotated_A(:,2) <= first_lamellae_upper;
    second_lamellae = rotated_A(:,2) >= second_lamellae_lower & rotated_A(:,2) <= second_lamellae_upper;
    first_center  = mean(rotated_A(first_lamellae, 2));
    second_center = mean(rotated_A(second_lamellae, 2));
    optional_third_lamellae = rotated_A(:,2) >= optional_third_lower & rotated_A(:,2) <= optional_third_upper;
    third_center = mean(rotated_A(optional_third_lamellae, 2));
    optional_fourth_lamellae = rotated_A(:,2) >= optional_fourth_lower & rotated_A(:,2) <= optional_fourth_upper;
    fourth_center = mean(rotated_A(optional_fourth_lamellae, 2));
    domain_spacing = (abs(second_center - first_center) + abs(third_center-second_center) + abs(fourth_center - third_center))/3;
    plane_SD=mean([std(rotated_A(first_lamellae, 2)),std(rotated_A(second_lamellae, 2)), std(rotated_A(optional_third_lamellae,2)),std(rotated_A(optional_fourth_lamellae,2))]);
else
    domain_spacing = abs(second_center - first_center);
end
box_dim=box./3;

fprintf('BOX DS X %.4g: BOX DS Y %.4g: BOX DS Z %.4g\n', box_dim);
disp(sprintf("DOMAIN SPACING : %g",domain_spacing))
disp(sprintf("SD FROM MEAN : %g",plane_SD))
fprintf('Lx %g: Ly %g: Lz %g:',box);

figure()
scatter3(A_pos(:,1), A_pos(:,2), A_pos(:,3), 'filled','blue')
xlabel('x')
ylabel('y')
zlabel('z')
ax = gca;                % get current axes
ax.LineWidth = 2;        % makes axis lines thicker (bolder)
ax.FontWeight = 'bold';  % makes tick labels bold
ax.FontSize = 14;

figure()
scatter3(rotated_A(:,1),rotated_A(:,2), rotated_A(:,3), 'filled','blue')
xlabel('x')
ylabel('y')
zlabel('z')
ax = gca;                % get current axes
ax.LineWidth = 2;        % makes axis lines thicker (bolder)
ax.FontWeight = 'bold';  % makes tick labels bold
ax.FontSize = 14;
legend(sprintf('Domain spacing = %.2f, SD = %.2f', domain_spacing, plane_SD), ...
       'Location', 'best');
second_lamellae_total=rotated_A(second_lamellae~=0,:);
max_X_plane = max(second_lamellae_total(:,1));
max_Y_plane = max(second_lamellae_total(:,2));
max_Z_plane = max(second_lamellae_total(:,3));
min_X_plane = min(second_lamellae_total(:,1));
min_Y_plane = min(second_lamellae_total(:,2));
min_Z_plane = min(second_lamellae_total(:,3));
lamellar_height = max_Z_plane - min_Z_plane;
lamellar_length = max_X_plane - min_X_plane;
X_dimension = [-box(1)/2 box(1)/2]+2;
Y_dimension = [-box(2)/2 box(2)/2];
Z_dimension = [-box(3)/2 box(3)/2];
second_lamellae_trimmed_X = ...
    second_lamellae_total(:,3) < Z_dimension(2) & ...
    second_lamellae_total(:,3) > Z_dimension(1) & ...
    second_lamellae_total(:,1) > X_dimension(1) & ...
    second_lamellae_total(:,1) < X_dimension(2);
saved_points = second_lamellae_total(second_lamellae_trimmed_X~=0,:);
mean_Y = mean(saved_points(:,2));
added_point_1 = [23 mean_Y -16];
added_point_2 = [23 mean_Y -7];
saved_points_middle = [saved_points; added_point_1; added_point_2]; 
figure()
scatter(saved_points(:,1),saved_points(:,3))
yline(-box(2)/2,'r')
yline(box(2)/2,'r')
xline(-box(1)/2,'r')
xlim([-box(1)/2 box(1)/2])
ylim([-box(2)/2 box(2)/2])

% figure()
% scatter(rotated_C(:,2), rotated_C(:,3), 'filled','blue')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% ax = gca;                % get current axes
% ax.LineWidth = 2;        % makes axis lines thicker (bolder)
% ax.FontWeight = 'bold';  % makes tick labels bold
% ax.FontSize = 14;

