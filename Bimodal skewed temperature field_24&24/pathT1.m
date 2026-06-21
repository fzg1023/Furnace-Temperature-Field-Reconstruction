clc;
close all;
clear;

% spatial extent
x = linspace(0, 24, 1000);
y = linspace(0, 24, 1000);
[X, Y] = meshgrid(x, y);

% Temperature model formula
T = 800 + 420 ./ (0.018 * ((X - 18).^2 + (Y - 15).^2) + 1)+ 360 ./ (0.015 * ((X - 6).^2 + (Y - 9).^2) + 1);


figure;
surf(X, Y, T, 'EdgeColor', 'none');
colorbar;
title('Temperature field model');
xlabel('X');
ylabel('Y');
zlabel('temperature');
colormap(jet); 

% Define grid division points
grid_x = [0, 6, 12, 18, 24]; % X
grid_y = [0, 6, 12, 18, 24]; % Y
% Initialize the grid average temperature matrix
grid_avg_temps = zeros(length(grid_x)-1, length(grid_y)-1);

% Calculate the average temperature of each grid
for i = 1:length(grid_x)-1
    for j = 1:length(grid_y)-1
       
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

disp('Grid average temperature matrix：');
disp(grid_avg_temps);

figure;
hold on;
title('Grid division and path display of temperature field');
xlabel('X');
ylabel('Y');
zlabel('temperature');
colormap(jet); 
format short;
% Define the position of the acoustic transducer
sensor_positions = [
    8, 0;   
    16, 0;  
    24, 8;  
    24, 16; 
    16, 24; 
    8, 24;  
    0, 16; 
    0, 8;   
];

% Draw the position of the acoustic transducer
scatter(sensor_positions(:, 1), sensor_positions(:, 2), 100, 'r', 'filled');
text(sensor_positions(:, 1), sensor_positions(:, 2), {'1', '2', '3', '4', '5', '6', '7', '8'}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'w', 'FontSize', 12);

% Define the path of the acoustic transducer
sensor_pairs = [
    1, 3; 
    1, 4; 
    1, 5; 
    1, 6; 
    1, 7; 
    1, 8; 
    2, 3; 
    2, 4; 
    2, 5; 
    2, 6; 
    2, 7; 
    2, 8; 
    3, 5; 
    3, 6; 
    3, 7; 
    3, 8; 
    4, 5; 
    4, 6; 
    4, 7; 
    4, 8; 
    5, 7; 
    5, 8; 
    6, 7; 
    6, 8; 
];

% Initialize intercept matrix
A = zeros(size(sensor_pairs, 1), numel(grid_avg_temps)); 

% Calculate the intercept length of each path for each grid
for k = 1:size(sensor_pairs, 1)
   
    start_sensor = sensor_positions(sensor_pairs(k, 1), :);
    end_sensor = sensor_positions(sensor_pairs(k, 2), :);
   
    for i = 1:length(grid_x)-1
        for j = 1:length(grid_y)-1
   
            x_start = grid_x(i);
            x_end = grid_x(i+1);
            y_start = grid_y(j);
            y_end = grid_y(j+1);
            
           
            l = RayTracing(start_sensor, end_sensor, [x_start, x_end, y_start, y_end]);
            
   
            A(k, (i-1)*(length(grid_y)-1) + j) = l;
        end
    end
end

disp('Intercept information matrix：');

row_labels = arrayfun(@(k) sprintf('path %d', k), 1:size(sensor_pairs, 1), 'UniformOutput', false);

col_labels = arrayfun(@(i, j) sprintf('grid (line%d,row%d)', i, j), ...
    repmat(1:length(grid_x)-1, 1, length(grid_y)-1), ...
    repelem(1:length(grid_y)-1, 1, length(grid_x)-1), 'UniformOutput', false);

A_table = array2table(A, 'RowNames', row_labels, 'VariableNames', col_labels);
disp(A_table);

% Path drawing
for i = 1:size(sensor_pairs, 1)
    start_sensor = sensor_positions(sensor_pairs(i, 1), :);
    end_sensor = sensor_positions(sensor_pairs(i, 2), :);
    plot([start_sensor(1), end_sensor(1)], [start_sensor(2), end_sensor(2)], 'r--', 'LineWidth', 1.5);
end

for i = 1:length(grid_x)
    plot([grid_x(i), grid_x(i)], [0, 24], 'k-', 'LineWidth', 1.5);
    plot([0, 24], [grid_y(i), grid_y(i)], 'k-', 'LineWidth', 1.5);
end


hold off;

x = grid_avg_temps(:);

% Calculate path integral temperature
b = A * x;

% Calculate the average temperature of the path (divided by the length of the path)
path_lengths = vecnorm(sensor_positions(sensor_pairs(:,2),:) - sensor_positions(sensor_pairs(:,1),:), 2, 2);
path_avg_temps = b ./ path_lengths;

%% visualization
format bank; 
disp('Average temperature of the path (℃):');
path_table = table(path_avg_temps, 'RowNames', arrayfun(@(k)sprintf('path %02d',k),1:size(sensor_pairs,1),'UniformOutput',false));
disp(path_table);

%% Add Error
noise_level = 0.03;

noise = noise_level .* path_avg_temps .* randn(size(path_avg_temps));

noisy_path_avg_temps = path_avg_temps + noise;

disp('Average temperature of the path after adding 3% error (℃):');
noisy_path_table = table(noisy_path_avg_temps, 'RowNames', arrayfun(@(k)sprintf('path %02d',k),1:size(sensor_pairs,1),'UniformOutput',false));
disp(noisy_path_table);

% Generate integrated temperature data with noise (for reconstruction)
b_noisy = noisy_path_avg_temps .* path_lengths;

%% Save data for reconstruction
save('TemperatureReconstructionData.mat', ...
     'A','grid_avg_temps','b','b_noisy','path_lengths', ...
     'grid_x','grid_y','sensor_positions','sensor_pairs','-v7');

save('RayGeometryData.mat', ...
    'grid_x', 'grid_y', ...
    'sensor_positions', 'sensor_pairs', ...
    'A', 'path_lengths', 'b', 'b_noisy', ...
    '-v7');





%% Ray tracing function (accurate calculation of intercept length)
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
    
    dir_x = sign(dx);
    dir_y = sign(dy);
    
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

