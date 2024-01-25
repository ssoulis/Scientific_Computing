% Parameters for the Poisson equation
n = 30; % Number of interior nodes for each value of y
m = 2;  % Number of boundary conditions
h = 1 / (n + 1); % Distance of successive nodes

% Define the grid
[X, Y] = meshgrid(linspace(0, 1, n + 2), linspace(0, 0.5, m + 2));
X = X(2:end-1, 2:end-1);
Y = Y(2:end-1, 2:end-1);

% Construct the Poisson matrix A and the right-hand side vector b
A = constructPoissonMatrix_1070263(X, Y, h);
b = A * ones(n * m, 1);

% Exact solution for Poisson equation
xsol = solveExactPoisson(X, Y);

% Set an upper bound on the number of iterations
maxiter = 4 * n;

% Test without preconditioning
[x_noprec, flag_noprec, relres_noprec, iter_noprec, resvec_noprec, errvec_noprec] = pcg_1070263(A, b, 1e-6, maxiter, 'none', [], xsol, 'resvec', [], 'errvec', []);

% Test with incomplete Cholesky preconditioning
[x_ic, flag_ic, relres_ic, iter_ic, resvec_ic, errvec_ic] = pcg_1070263(A, b, 1e-6, maxiter, 'ichol', [], xsol, 'resvec', [], 'errvec', []);

% Test with custom preconditioning
[x_custom, flag_custom, relres_custom, iter_custom, resvec_custom, errvec_custom] = pcg_1070263(A, b, 1e-6, maxiter, 'custom', [], xsol, 'resvec', [], 'errvec', []);

% Plotting
figure;

subplot(2, 1, 1);
plot(0:iter_noprec, resvec_noprec / norm(b), '-o', 'DisplayName', 'No Preconditioning');
hold on;
plot(0:iter_ic, resvec_ic / norm(b), '-o', 'DisplayName', 'IC(0) Preconditioning');
plot(0:iter_custom, resvec_custom / norm(b), '-o', 'DisplayName', 'Custom Preconditioning');
title('Relative Residual');
xlabel('Iteration');
ylabel('Relative Residual');
legend('Location', 'Best');
grid on;

subplot(2, 1, 2);
plot(0:iter_noprec, errvec_noprec, '-o', 'DisplayName', 'No Preconditioning');
hold on;
plot(0:iter_ic, errvec_ic, '-o', 'DisplayName', 'IC(0) Preconditioning');
plot(0:iter_custom, errvec_custom, '-o', 'DisplayName', 'Custom Preconditioning');
title('Relative Error');
xlabel('Iteration');
ylabel('Relative Error');
legend('Location', 'Best');
grid on;

% Display the plots
drawnow;
