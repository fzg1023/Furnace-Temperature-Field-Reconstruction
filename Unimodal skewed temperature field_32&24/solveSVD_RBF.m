function [x_hat, w_hat] = solveSVD_RBF(A, b, Phi)

    A_rbf    = A * Phi;            % 24×N
    w_hat    = pinv(A_rbf) * b;    % SVD 
    x_hat    = Phi * w_hat;        % 16×1
end
