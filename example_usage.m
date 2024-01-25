% Example usage:
D = diag([2, 3, 4]);
L = tril(ones(3) - eye(3), -1);
U = triu(ones(3) - eye(3), 1);

A = D - L - U;
B = rand(3, 1);

% Call the function from the other script
x = solve_triangular_system(A, B);

% Display the result
disp(x);