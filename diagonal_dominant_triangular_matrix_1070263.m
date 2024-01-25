% Example usage with a diagonal-dominant triangular matrix
D = diag([10, 20, 30]);
L = tril(ones(3) - eye(3), -1);
U = triu(ones(3) - eye(3), 1);

A = D - L - U;
B = rand(3, 1);

% Set a small epsilon value for regularization
epsilon = 1e-6;

% Call the experiment function with defined epsilon
solve_triangular_system_experiment(A, B, epsilon);

% Additional experiments with diagonal-dominant triangular matrices

% Experiment 1: 5x5 matrix
A1 = diag(1:5) + tril(ones(5,5), -1);
B1 = rand(5, 1);
epsilon1 = 1e-6;
solve_triangular_system_experiment(A1, B1, epsilon1);

% Experiment 2: 7x7 matrix
A2 = diag(1:7) + tril(ones(7,7), -1);
B2 = rand(7, 1);
epsilon2 = 1e-6;
solve_triangular_system_experiment(A2, B2, epsilon2);

% Experiment 3: 10x10 matrix
A3 = diag(1:10) + tril(ones(10,10), -1);
B3 = rand(10, 1);
epsilon3 = 1e-6;
solve_triangular_system_experiment(A3, B3, epsilon3);