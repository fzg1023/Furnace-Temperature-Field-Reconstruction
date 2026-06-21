%% 32×24 
%  X：32m，Y：24m
clc;
close all;
clear;

x = linspace(0, 32, 1000);   % X
y = linspace(0, 24, 1000);   % Y
[X, Y] = meshgrid(x, y);

% Temperature model formula
T = 800 + 600 ./ (0.018 * ((X -20).^2 + (Y - 15).^2) + 1);


figure('Position', [100, 100, 800, 600]);  
imagesc(x, y, T);
axis xy equal tight;          
colorbar;
title('Temperature field model (32m × 24m)', 'FontSize', 14);
xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
colormap(jet);


xlim([0 32]);
ylim([0 24]);
set(gca, 'FontSize', 12);
set(gca, 'DataAspectRatio', [1 1 1]);  


%% Define grid division points
grid_x = [0, 8, 16, 24, 32];     
grid_y = [0, 6, 12, 18, 24];     

% Initialize the grid average temperature matrix
n_grids_x = length(grid_x) - 1;  
n_grids_y = length(grid_y) - 1; 
grid_avg_temps = zeros(n_grids_x, n_grids_y);

% Calculate the average temperature of each grid
for i = 1:n_grids_x
    for j = 1:n_grids_y
        
        x_start = grid_x(i);
        x_end = grid_x(i+1);
        y_start = grid_y(j);
        y_end = grid_y(j+1);
        

        mask = (X >= x_start) & (X < x_end) & (Y >= y_start) & (Y < y_end);
        grid_temps = T(mask);
        
     
        grid_avg_temps(i, j) = mean(grid_temps(:));
    end
end

format bank;
disp('Grid average temperature matrix（4×4）：');
disp(grid_avg_temps);


fprintf('\nmesh generation：X=%s，Y=%s\n', ...
        mat2str(grid_x), mat2str(grid_y));
fprintf('grid number：%d × %d = %d\n', n_grids_x, n_grids_y, numel(grid_avg_temps));


sensor_positions = [
     8,  0;   
    16,  0;   
    24,  0;   
    32,  8;  
    32, 16;   
    24, 24;   
    16, 24;   
     8, 24;   
     0, 16;   
     0,  8;   
];

%% 4. Define sound wave path
sensor_pairs = [];
for i = 1:size(sensor_positions, 1)
    for j = i+1:size(sensor_positions, 1)
        
        xi = sensor_positions(i, 1);
        yi = sensor_positions(i, 2);
        xj = sensor_positions(j, 1);
        yj = sensor_positions(j, 2);
        
        
        
        if yi == 0 && yj == 0
            continue;
        end
       
        if yi == 24 && yj == 24
            continue;
        end
       
        if xi == 0 && xj == 0
            continue;
        end
        
        if xi == 32 && xj == 32
            continue;
        end
        
        
        sensor_pairs = [sensor_pairs; i, j];
    end
end

n_paths = size(sensor_pairs, 1);
n_sensors = size(sensor_positions, 1);
fprintf('Number of sensors：%d\n', n_sensors);
fprintf('Theoretical maximum number of paths：C(%d,2) = %d\n', n_sensors, nchoosek(n_sensors, 2));
fprintf('Actual number of effective paths：%d\n', n_paths);

%%  Draw sensor position and path
figure('Position', [100, 100, 800, 600]);
hold on;
title('Grid division of temperature field and acoustic path (32m × 24m)', 'FontSize', 14);
xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
colormap(jet);


for i = 1:length(grid_x)
    plot([grid_x(i), grid_x(i)], [0, 24], 'k-', 'LineWidth', 1.5);
end
for j = 1:length(grid_y)
    plot([0, 32], [grid_y(j), grid_y(j)], 'k-', 'LineWidth', 1.5);
end


for i = 1:size(sensor_pairs, 1)
    start_sensor = sensor_positions(sensor_pairs(i, 1), :);
    end_sensor = sensor_positions(sensor_pairs(i, 2), :);
    plot([start_sensor(1), end_sensor(1)], [start_sensor(2), end_sensor(2)], ...
         'r--', 'LineWidth', 1.2);
end

% Draw sensor position
scatter(sensor_positions(:, 1), sensor_positions(:, 2), 150, 'r', 'filled');
labels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
for i = 1:n_sensors
    x_shift = sensor_positions(i, 1);
    y_shift = sensor_positions(i, 2);
   
    offset = 1.2;
    if y_shift == 0
        y_shift = y_shift - offset;
    elseif y_shift == 24
        y_shift = y_shift + offset;
    elseif x_shift == 0
        x_shift = x_shift - offset;
    elseif x_shift == 32
        x_shift = x_shift + offset;
    end
    text(x_shift, y_shift, labels{i}, ...
        'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

axis equal tight;
xlim([0 32]);
ylim([0 24]);
xticks(grid_x);
yticks(grid_y);
set(gca, 'FontSize', 12);
box on;
hold off;

%%  Calculate intercept matrix A
n_grids = numel(grid_avg_temps);  % 16个网格
A = zeros(n_paths, n_grids);

fprintf('\nCalculating intercept matrix...\n');
tic;
%% Construction of Revised A Matrix
for k = 1:n_paths
    start_sensor = sensor_positions(sensor_pairs(k, 1), :);
    end_sensor = sensor_positions(sensor_pairs(k, 2), :);
    
    for i = 1:n_grids_x
        for j = 1:n_grids_y
            x_start = grid_x(i);
            x_end = grid_x(i+1);
            y_start = grid_y(j);
            y_end = grid_y(j+1);
            
            l = RayTracing(start_sensor, end_sensor, [x_start, x_end, y_start, y_end]);
            
           
            col = (j-1) * n_grids_x + i; 
            A(k, col) = l;
        end
    end
end
elapsed = toc;
fprintf('The intercept matrix calculation is complete! time-consuming %.2f s\n', elapsed);



fprintf('  A: %d × %d (paths × grids)\n', size(A,1), size(A,2));


nz_ratio = nnz(A) / numel(A) * 100;
fprintf('  Non zero element ratio：%.2f%%\n', nz_ratio);

% 检查每列（每个网格）被穿过的路径数
grid_coverage = sum(A > 0, 1);
fprintf('  The number of paths traversed by each grid：%s\n', mat2str(grid_coverage));

% 显示前几条路径的截距
disp('The intercept information of the first 5 paths：');
for k = 1:min(5, n_paths)
    fprintf('path %d (sensor%d→%d): ', k, sensor_pairs(k,1), sensor_pairs(k,2));
    for col = 1:n_grids
        if A(k, col) > 0
            fprintf('grid%d=%.2f  ', col, A(k, col));
        end
    end
    fprintf('\n');
end

%% Calculate path integral temperature and average temperature
x_true = grid_avg_temps(:);


b = A * x_true;

% Calculate path length
path_lengths = zeros(n_paths, 1);
for k = 1:n_paths
    start_pt = sensor_positions(sensor_pairs(k, 1), :);
    end_pt = sensor_positions(sensor_pairs(k, 2), :);
    path_lengths(k) = norm(end_pt - start_pt);
end


path_avg_temps = b ./ path_lengths;


format bank;
disp('Average temperature of the path (top 10):');
for k = 1:min(10, n_paths)
    fprintf('  path %02d (sensor%d→%d): length=%.2fm, average temperature=%.2f K\n', ...
            k, sensor_pairs(k,1), sensor_pairs(k,2), ...
            path_lengths(k), path_avg_temps(k));
end

%% 1 Add Error
noise_level = 0.03;  


noise = noise_level .* path_avg_temps .* randn(size(path_avg_temps));
noisy_path_avg_temps = path_avg_temps + noise;



for k = 1:min(10, n_paths)
    fprintf('  path %02d: true=%.2f K, Containing noise=%.2f K, noise=%.2f K\n', ...
            k, path_avg_temps(k), noisy_path_avg_temps(k), ...
            noisy_path_avg_temps(k)-path_avg_temps(k));
end


b_noisy = noisy_path_avg_temps .* path_lengths;

%% 11. save data
save('TemperatureReconstructionData_32x24.mat', ...
     'A', 'grid_avg_temps', 'b', 'b_noisy', 'path_lengths', ...
     'grid_x', 'grid_y', 'sensor_positions', 'sensor_pairs', '-v7');

save('RayGeometryData_32x24.mat', ...
    'grid_x', 'grid_y', ...
    'sensor_positions', 'sensor_pairs', ...
    'A', 'path_lengths', 'b', 'b_noisy', ...
    '-v7');

fprintf('\nThe data has been saved to:\n');
fprintf('  TemperatureReconstructionData_32x24.mat\n');
fprintf('  RayGeometryData_32x24.mat\n');

%%  total
fprintf('\n═══════════════════════════════════\n');
fprintf('  Summary of Boiler Structural Parameters\n');
fprintf('═══════════════════════════════════\n');
fprintf('  Boiler size：%.0fm × %.0fm (X×Y)\n', grid_x(end), grid_y(end));
fprintf('  mesh generation：%d × %d = %d grids\n', ...
        n_grids_x, n_grids_y, n_grids);
fprintf('  grid size：%.0fm × %.0fm\n', ...
        grid_x(2)-grid_x(1), grid_y(2)-grid_y(1));
fprintf('  Number of sensors：%d\n', n_sensors);
fprintf('  Number of sound wave paths：%d\n', n_paths);
fprintf('  Intercept matrix A：%d × %d\n', size(A, 1), size(A, 2));
fprintf('  Non zero element ratio：%.2f%%\n', nz_ratio);
fprintf('  Error level：%.0f%%\n', noise_level * 100);
fprintf('═══════════════════════════════════\n');

%% Intercept calculation
function l = RayTracing(start_sensor, end_sensor, grid_bounds)
   
    x_min = grid_bounds(1);
    x_max = grid_bounds(2);
    y_min = grid_bounds(3);
    y_max = grid_bounds(4);
    

    x1 = start_sensor(1);
    y1 = start_sensor(2);
    x2 = end_sensor(1);
    y2 = end_sensor(2);
    dx = x2 - x1;
    dy = y2 - y1;
    

    t_min = 0;
    t_max = 1;
    

    if dx ~= 0
        t_xmin = (x_min - x1) / dx;
        t_xmax = (x_max - x1) / dx;
        t_x1 = min(t_xmin, t_xmax);
        t_x2 = max(t_xmin, t_xmax);
        t_min = max(t_min, t_x1);
        t_max = min(t_max, t_x2);
    else
        if x1 < x_min || x1 >= x_max
            l = 0;
            return;
        end
    end
    

    if dy ~= 0
        t_ymin = (y_min - y1) / dy;
        t_ymax = (y_max - y1) / dy;
        t_y1 = min(t_ymin, t_ymax);
        t_y2 = max(t_ymin, t_ymax);
        t_min = max(t_min, t_y1);
        t_max = min(t_max, t_y2);
    else
        if y1 < y_min || y1 >= y_max
            l = 0;
            return;
        end
    end
    
   
    if t_min < t_max && t_max > 0 && t_min < 1
       
        t_enter = max(t_min, 0);
        t_exit = min(t_max, 1);
        l = norm([dx, dy]) * (t_exit - t_enter);
    else
        l = 0;
    end
end