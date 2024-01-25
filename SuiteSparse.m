A = Problem;

% Define right-hand side vector
b = rand(size(A, 1), 1);

% Set preconditioning options
preconditioner = 'ichol';  % Change to 'none', 'custom', or other options as needed

% Set an upper bound on the number of iterations
maxiter = 4 * size(A, 1);

% Run PCG with modified function
[x, flag, relres, iter, resvec, errvec] = pcg_myid_1070263(A, b, 1e-6, maxiter, preconditioner, [], [], 'resvec', [], 'errvec', []);

% Plotting
figure;
subplot(2, 1, 1);
plot(0:iter, resvec / norm(b), '-o', 'DisplayName', 'Relative Residual');
title('Relative Residual');
xlabel('Iteration');
ylabel('Relative Residual');
legend('Location', 'Best');
grid on;

subplot(2, 1, 2);
plot(0:iter, errvec, '-o', 'DisplayName', 'Relative Error');
title('Relative Error');
xlabel('Iteration');
ylabel('Relative Error');
legend('Location', 'Best');
grid on;

% Display the plots
drawnow;
