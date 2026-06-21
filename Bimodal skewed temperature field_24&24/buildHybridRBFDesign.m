function Phi = buildHybridRBFDesign(cx, cy, sigma1, sigma2, usePR)
%Constructing a mixed kernel function (Gaussian+inverse multiple quadratic) design matrix
%   cx, cy   ：Gaussian kernel center (length of 16)
%   sigma1   ：Scale parameters of Gaussian kernel
%   sigma2   ：Scale parameters of inverse multiple quadratic kernels
%   usePR    ：Enable polynomial regeneration term
%   Phi      ：Design Matrix（ 16 × (16+4+3) = 16 × 23）

    n_g = numel(cx);    %  16
    m = n_g;            
    X = cx;
    Y = cy;

    % 1. Gaussia（16 ）
    Phi_g = zeros(m, n_g);
    for i = 1:n_g
        r2 = (X - cx(i)).^2 + (Y - cy(i)).^2;
        Phi_g(:, i) = exp(-r2 / (2 * sigma1^2));
    end

    % 2. inverse multiple quadratic（4 ）
    imq_centers = [6 6; 6 18; 18 6; 18 18]; 
    n_f = size(imq_centers, 1);
    Phi_f = zeros(m, n_f);
    for i = 1:n_f
        xc = imq_centers(i,1);
        yc = imq_centers(i,2);
        r2 = (X - xc).^2 + (Y - yc).^2;
        Phi_f(:, i) = 1 ./ sqrt(r2 + sigma2^2);  % IMQ
    end

    % 3. polynomial regeneration（3 ）
    if usePR
        Phi_pr = [ones(m,1), X, Y];  
        Phi = [Phi_g, Phi_f, Phi_pr];  
    else
        Phi = [Phi_g, Phi_f]; 
    end
end
