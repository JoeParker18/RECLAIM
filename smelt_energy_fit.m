% This code is meant strictly just to solve the matrix to get E(0) for more
% accurate calculations for approximating the time for SMELT

%% A and B matrices

A = [1, 19.430, 1/2;
     1, 24.490, 1/2;
     1, 51.810, 1/2;
     1, 18.330, 1/2;
     1, 92.160, 1/2.5;
     1, 151.76, 1/2.5;
     1, 163.40, 1/2.5;
     1, 247.59, 1/2.5;
     1, 161.60, 1/3;
     1, 430.46, 1/3;
     1, 104.72, 1/3];

b = [480;
     420;
     480;
     552;
     600;
     900;
     900;
     1350;
     720;
     1800;
     798;];

%% Solving the matrix
coef = A\b;
E_0 = coef(1); e = coef(2); k = coef(3);

% --- solve ---
fprintf('rows: %d, rank = %d (need 3)\n', size(A,1), rank(A));
coef = A \ b;
E0 = coef(1);  e = coef(2);  k = coef(3);

% --- fit quality ---
pred = A*coef;
R2 = 1 - sum((b-pred).^2)/sum((b-mean(b)).^2);
cc = corrcoef(A(:,2), A(:,3));

fprintf('E_run = %.1f + %.3f*M + %.1f/P\n', E0, e, k);
fprintf('R2 = %.4f\n', R2);
fprintf('corr(M, 1/P) = %.3f\n', cc(1,2));