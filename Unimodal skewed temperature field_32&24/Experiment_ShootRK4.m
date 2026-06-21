%% Experiment_ShootRK4_32x24.m
% ------------------------------------------------------------------
clc; clear; close all;

%% data
load TemperatureReconstructionData_32x24.mat       
b_obs = b_noisy;                                   

%%  Sound velocity&prior temperature function
Z    = 19.63;                                      % c = Z*sqrt(T)
Tfun = @(x,y) 800 + 360 ./ (0.018*((x-12).^2 + (y-15).^2) + 1);


%%  A_curve 
[A_curve, PathNodes] = buildAcurve_shootRK4( ...
        Tfun, grid_x, grid_y, sensor_positions, sensor_pairs);

A = A_curve;                 

fprintf('Curved ray geometry matrix：%d paths，%d grids\n', size(A,1), size(A,2));

%% shape parameter 
[cx, cy] = ndgrid( (grid_x(1:end-1)+grid_x(2:end))/2 , ...
                   (grid_y(1:end-1)+grid_y(2:end))/2 );
cx = cx(:); cy = cy(:); usePR = true;

D = pdist([cx cy]); d_avg = mean(D);
sig1 = linspace(d_avg/2, 2*d_avg, 50); rmse1 = zeros(size(sig1));

[xQ, yQ] = meshgrid(linspace(grid_x(1), grid_x(end), 134), ...
                    linspace(grid_y(1), grid_y(end), 100));
Tq_true = 800 + 600 ./ (0.018*((xQ-20).^2 + (yQ-15).^2) + 1);


for k = 1:numel(sig1)
    Phi_k = buildRBFDesign(cx, cy, sig1(k), usePR);
    [~, w_k] = solveSVD_RBF(A, b_obs, Phi_k);
    Tq_est = zeros(size(xQ));
    for p = 1:numel(cx)
        G = exp(-((xQ - cx(p)).^2 + (yQ - cy(p)).^2)/(2*sig1(k)^2));
        Tq_est = Tq_est + w_k(p)*G;
    end
    Tq_est = Tq_est + w_k(end-2) + w_k(end-1)*xQ + w_k(end)*yQ;
    rmse1(k) = sqrt(mean((Tq_est(:) - Tq_true(:)).^2));
end
sigma_c =0.98 * sig1(rmse1 == min(rmse1));
sigma_f = 1.05 * sigma_c;
fprintf('σ_c = %.3f,  σ_f = %.3f\n', sigma_c, sigma_f);

%% Design Matrix（16 + 4 + 3PR = 23）
Phi = buildHybridRBFDesign(cx, cy, sigma_c, sigma_f, usePR);

%%  Wall temperature constraint (10 sensor positions)
wall_points = sensor_positions;  
Nw = size(wall_points, 1);


imq_centers = [8 6; 24 6; 8 18; 24 18];

L_c = zeros(Nw, numel(cx));
L_f = zeros(Nw, 4);

for i = 1:Nw
    xw = wall_points(i,1); yw = wall_points(i,2);

    for j = 1:numel(cx)
        L_c(i,j) = exp(-((xw-cx(j))^2 + (yw-cy(j))^2)/(2*sigma_c^2));
    end
  
    for j = 1:4
        r2 = (xw-imq_centers(j,1)).^2 + (yw-imq_centers(j,2)).^2;
        L_f(i,j) = 1/sqrt(r2 + sigma_f^2);
    end
end

L_PR   = [ones(Nw,1) wall_points];
L_wall = [L_c L_f L_PR];


T_wall = 800 + 600 ./ (0.018*((wall_points(:,1)-20).^2 + (wall_points(:,2)-15).^2) + 1);


fprintf('Wall temperature constraint point：%d（sensor position）\n', Nw);

%% L-curve -- α
alphas = logspace(-3, 3, 30); 
data_f = zeros(size(alphas));
prior_f = data_f; 
rank_keep = data_f;

for k = 1:length(alphas)
    a = alphas(k);
    M = [A*Phi; sqrt(a)*L_wall];
    b = [b_obs; sqrt(a)*T_wall];
    [U, S, V] = svd(M, 'econ'); 
    s = diag(S); 
    tol = max(s)*3e-4;
    r = find(s > tol, 1, 'last');
    w = V(:,1:r) * ((U(:,1:r)' * b) ./ s(1:r));
    data_f(k)  = norm(A*Phi*w - b_obs);
    prior_f(k) = norm(L_wall*w - T_wall);
    rank_keep(k) = r;
end

d = (log10(data_f)-mean(log10(data_f))).^2 + ...
    (log10(prior_f)-mean(log10(prior_f))).^2;
alpha_opt = alphas(d==min(d)) + 1000;
fprintf('α_opt = %.3e  (rank = %d)\n', alpha_opt, rank_keep(d==min(d)));

%% 7. Request weight
M_aug = [A*Phi; sqrt(alpha_opt)*L_wall];
b_aug = [b_obs; sqrt(alpha_opt)*T_wall];
[U, S, V] = svd(M_aug, 'econ'); 
s = diag(S); 
tol = max(s)*3e-4;
r = find(s > tol, 1, 'last');
w_hat = V(:,1:r) * ((U(:,1:r)' * b_aug) ./ s(1:r));

fprintf('Final TSVD truncation order：r=%d\n', r);

%% rebuild
ngr_x = 500;
ngr_y = 375;
[xFine, yFine] = meshgrid(linspace(grid_x(1), grid_x(end), ngr_x), ...
                          linspace(grid_y(1), grid_y(end), ngr_y));
T_est = zeros(size(xFine));


for p = 1:numel(cx)
    Gc = exp(-((xFine-cx(p)).^2 + (yFine-cy(p)).^2)/(2*sigma_c^2));
    T_est = T_est + w_hat(p)*Gc;
end


for j = 1:4
    cx_imq = imq_centers(j,1); 
    cy_imq = imq_centers(j,2);
    Gf = 1 ./ sqrt((xFine-cx_imq).^2 + (yFine-cy_imq).^2 + sigma_f^2);
    T_est = T_est + w_hat(numel(cx)+j)*Gf;
end


offset = numel(cx) + 4;
T_est  = T_est + w_hat(offset+1) + w_hat(offset+2)*xFine + w_hat(offset+3)*yFine;

%% 9. error
T_true = 800 + 600 ./ (0.018*((xFine-20).^2 + (yFine-15).^2) + 1) ;


absErr = abs(T_est - T_true);
MaxAE  = max(absErr(:));
MRE    = mean(absErr(:) ./ T_true(:));
RMSRE  = sqrt(mean((absErr(:) ./ T_true(:)).^2));

%% Boundary zone and center zone error
wth_x = grid_x(2) - grid_x(1);    % X
wth_y = grid_y(2) - grid_y(1);    % Y 

isBd = (xFine <= grid_x(1) + wth_x) | ...
       (xFine >= grid_x(end) - wth_x) | ...
       (yFine <= grid_y(1) + wth_y) | ...
       (yFine >= grid_y(end) - wth_y);

% Eb
bdAbs = absErr(isBd); 
bdRel = bdAbs ./ T_true(isBd);
MaxAE_bd = max(bdAbs); 
MRE_bd = mean(bdRel);
RMSE_bd = sqrt(mean(bdAbs.^2)); 
CV_RMSE_bd = RMSE_bd / mean(T_true(isBd)) * 100;

% Ea
centerAbs = absErr(~isBd);
centerRel = centerAbs ./ T_true(~isBd);
MaxAE_center = max(centerAbs);
MRE_center = mean(centerRel);
RMSE_center = sqrt(mean(centerAbs.^2));
CV_RMSE_center = RMSE_center / mean(T_true(~isBd)) * 100;
% ---------- 5 true & error -----------------------------------------------
T_true = 800 + 600 ./ ( 0.018*((xFine-20).^2 + (yFine-15).^2) + 1 );


signedErr = T_est - T_true;           


[MaxErr, idx_max] = max(signedErr(:));    
[MinErr, idx_min] = min(signedErr(:));   

if abs(MaxErr) >= abs(MinErr)
    MaxSignedErr = MaxErr;
else
    MaxSignedErr = MinErr;
end

[x_max, y_max] = ind2sub(size(signedErr), idx_max);
[x_min, y_min] = ind2sub(size(signedErr), idx_min);

fprintf('Statistics of signed errors：\n');
fprintf('  Maximum positive error = %.3f K  (position: x=%.2f, y=%.2f)\n', MaxErr, xFine(idx_max), yFine(idx_max));
fprintf('  Maximum negative error = %.3f K  (position: x=%.2f, y=%.2f)\n', MinErr, xFine(idx_min), yFine(idx_min));
fprintf('  Global maximum error = %.3f K\n', MaxSignedErr);
%% 10. printf
fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║  error (32m×24m)       ║\n');
fprintf('╠══════════════════════════════════════════════════╣\n');
fprintf('║  σ_c=%.2f, σ_f=%.2f, α=%.2e, r=%d              ║\n', sigma_c, sigma_f, alpha_opt, r);
fprintf('╠══════════════════════════════════════════════════╣\n');
fprintf('║  【Global】                                ║\n');
fprintf('║    MaxAE  = %.3f K                              ║\n', MaxAE);
fprintf('║    MRE    = %.4f (%.2f%%)                        ║\n', MRE, 100*MRE);
fprintf('║    RMSRE  = %.4f (%.2f%%)                        ║\n', RMSRE, 100*RMSRE);
fprintf('╠══════════════════════════════════════════════════╣\n');
fprintf('║  【Boundary zone error】(X=%.0fm, Y=%.0fm)           ║\n', wth_x, wth_y);
fprintf('║    MaxAE  = %.3f K                              ║\n', MaxAE_bd);
fprintf('║    MRE    = %.4f (%.2f%%)                        ║\n', MRE_bd, 100*MRE_bd);
fprintf('║    CV(RMSE) = %.2f %%                            ║\n', CV_RMSE_bd);
fprintf('╠══════════════════════════════════════════════════╣\n');
fprintf('║  【Center error】                       ║\n');
fprintf('║    MaxAE  = %.3f K                              ║\n', MaxAE_center);
fprintf('║    MRE    = %.4f (%.2f%%)                        ║\n', MRE_center, 100*MRE_center);
fprintf('║    CV(RMSE) = %.2f %%                            ║\n', CV_RMSE_center);
fprintf('╚══════════════════════════════════════════════════╝\n');

%% 11. visualization
cmin = min(T_true(:)); 
cmax = max(T_true(:));

% Real Field
figure;
imagesc(linspace(0,32,ngr_x), linspace(0,24,ngr_y), T_true);
axis xy equal tight;
xlim([0 32]); ylim([0 24]);
xticks(grid_x); yticks(grid_y);
xlabel('x/m', 'FontSize', 16);
ylabel('y/m', 'FontSize', 16);
title('', 'FontSize', 16);
set(gca, 'FontSize', 16);
caxis([cmin cmax]); colorbar; colormap(jet);

% Reconstruction site
figure;
imagesc(linspace(0,32,ngr_x), linspace(0,24,ngr_y), T_est);
axis xy equal tight;
xlim([0 32]); ylim([0 24]);
xticks(grid_x); yticks(grid_y);
xlabel('x/m', 'FontSize', 16);
ylabel('y/m', 'FontSize', 16);
title(sprintf('', sigma_c, sigma_f), 'FontSize', 16);
set(gca, 'FontSize', 16);
caxis([cmin cmax]); colorbar; colormap(jet);

% error
figure;
signedErr = T_est - T_true;
surf(xFine, yFine, signedErr, 'EdgeColor', 'none');
axis tight;
xlim([0 32]); ylim([0 24]);
xticks(grid_x); yticks(grid_y);
xlabel('x/m', 'FontSize', 16);
ylabel('y/m', 'FontSize', 16);
zlabel('Error (K)', 'FontSize', 16);
title(sprintf('', sigma_c), 'FontSize', 16);
set(gca, 'FontSize', 16);
maxAbsErr = max(abs(signedErr(:)));
caxis([-maxAbsErr, maxAbsErr]);
colorbar; colormap(jet);
view(-30, 15);
shading interp; lighting gouraud;



%% 12. Schematic diagram of curved ray path
figure('Position', [100, 50, 900, 650]); 
hold on;
title('', 'FontSize', 14);
xlabel('x/m', 'FontSize', 12); 
ylabel('y/m', 'FontSize', 12); 
axis equal tight;


for i = 1:numel(grid_x)
    plot([grid_x(i) grid_x(i)], [grid_y(1) grid_y(end)], 'k:', 'LineWidth', 1);
end
for j = 1:numel(grid_y)
    plot([grid_x(1) grid_x(end)], [grid_y(j) grid_y(j)], 'k:', 'LineWidth', 1);
end


for k = 1:numel(PathNodes)
    plot(PathNodes{k}(:,1), PathNodes{k}(:,2), ...
         'Color', [1 0.6 0.6], 'LineWidth', 1.2);
end


scatter(sensor_positions(:,1), sensor_positions(:,2), 100, 'b', 'filled');


labels = {'A','B','C','D','E','F','G','H','I','J'};
offset = 1.2;
for i = 1:size(sensor_positions,1)
    x_shift = sensor_positions(i,1);
    y_shift = sensor_positions(i,2);
    if y_shift == 0
        y_shift = y_shift - offset-0.3;
    elseif y_shift == 24
        y_shift = y_shift + offset;
    elseif x_shift == 0
        x_shift = x_shift - offset;
    elseif x_shift == 32
        x_shift = x_shift + offset;
    end
    text(x_shift, y_shift, labels{i}, ...
        'Color', [0.12 0.60 1.00], 'FontSize', 14, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

xlim([0 32]); ylim([0 24]);
xticks(grid_x); yticks(grid_y);
set(gca, 'FontSize', 14);
box on; 
hold off;



%% (shoot+RK4)
function [A_curve, PathNodes] = buildAcurve_shootRK4( ...
                        Tfun, grid_x, grid_y, sensor_pos, sensor_pairs)
    nPath = size(sensor_pairs, 1);
    nGrid = (numel(grid_x)-1)*(numel(grid_y)-1);
    A_curve   = zeros(nPath, nGrid);
    PathNodes = cell(nPath, 1);
    
    fprintf('Calculating curved ray path...\n');
    for k = 1:nPath
        p1 = sensor_pos(sensor_pairs(k,1), :);
        p2 = sensor_pos(sensor_pairs(k,2), :);
        [nodes, segLen] = RayShootRK4(p1, p2, Tfun);
        PathNodes{k} = nodes;
 
        for s = 1:numel(segLen)
            mid = 0.5*(nodes(s,:) + nodes(s+1,:));
            i = find(grid_x <= mid(1), 1, 'last');
            j = find(grid_y <= mid(2), 1, 'last');
            if i>=1 && i<numel(grid_x) && j>=1 && j<numel(grid_y)
                col = (i-1)*(numel(grid_y)-1) + j;
                A_curve(k, col) = A_curve(k, col) + segLen(s);
            end
        end
        if mod(k, 10) == 0
            fprintf('  Processed %d/%d paths \n', k, nPath);
        end
    end
    fprintf('Curved ray path calculation completed!\n');
end

%% --- RayShootRK4 
function [XY, segLen] = RayShootRK4(p1, p2, Tfun)
    rot = false;
    if abs(p2(1)-p1(1)) < abs(p2(2)-p1(2))
        rot = true; 
        p1 = p1([2 1]); 
        p2 = p2([2 1]);
    end
    x1 = p1(1); x2 = p2(1);  
    y1 = p1(2); y2 = p2(2);
    L  = x2 - x1;
    N  = 800; 
    h = L/N;
    v1 = (y2-y1)/L;
    v2 = v1*1.1 + 1e-4;
    F1 = shoot(v1); 
    F2 = shoot(v2);
    iter = 0; 
    maxIt = 30; 
    tol = 1e-4;
    
    while abs(F2) > tol && iter < maxIt
        v3 = v2 - F2*(v2-v1)/(F2-F1+eps);
        v1 = v2; F1 = F2; 
        v2 = v3; F2 = shoot(v2);
        iter = iter + 1;
    end
    
    [~, y] = ivpSolve(v2);
    x = linspace(x1, x2, N+1).';
    if rot
        XY = [y x];
    else
        XY = [x y];
    end
    d = diff(XY, 1, 1); 
    segLen = sqrt(sum(d.^2, 2));


    function F = shoot(v0)
        [~, ytmp] = ivpSolve(v0);
        F = ytmp(end) - y2;
    end
    
    function [x, y] = ivpSolve(v0)
        x = linspace(x1, x2, N+1).';
        y = zeros(N+1, 1); 
        v = zeros(N+1, 1);
        y(1) = y1; 
        v(1) = v0;
        for i = 1:N
            xi = x(i); yi = y(i); vi = v(i);
            k1y = vi;
            k1v = rhs(xi, yi, vi);
            k2y = vi + 0.5*h*k1v;
            k2v = rhs(xi+0.5*h, yi+0.5*h*k1y, vi+0.5*h*k1v);
            k3y = vi + 0.5*h*k2v;
            k3v = rhs(xi+0.5*h, yi+0.5*h*k2y, vi+0.5*h*k2v);
            k4y = vi + h*k3v;
            k4v = rhs(xi+h, yi+h*k3y, vi+h*k3v);
            y(i+1) = yi + h*(k1y+2*k2y+2*k3y+k4y)/6;
            v(i+1) = vi + h*(k1v+2*k2v+2*k3v+k4v)/6;
        end
    end
    
    function f = rhs(x, y, v)
        hgrad = 1e-3;
        dTx = (Tfun(x+hgrad, y) - Tfun(x-hgrad, y))/(2*hgrad);
        dTy = (Tfun(x, y+hgrad) - Tfun(x, y-hgrad))/(2*hgrad);
        f = 0.5*(1+v^2)*(dTx - v*dTy)/Tfun(x, y);
    end
end