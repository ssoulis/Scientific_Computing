function experimentMatrixProperties(n)
    % Initialize the diagonal components with random values
    D = diag(randi([1, 5], n, 1)); % Main diagonal
    L = diag(randi([1, 5], n-1, 1), -1); % Lower diagonal
    U = diag(randi([1, 5], n-1, 1), 1); % Upper diagonal
    A = D - L - U;

    % Check if the matrix is positive definite
    try
        chol(A);
        disp('Matrix is positive definite.');
    catch
        disp('Matrix is not positive definite.');
    end

    % Initialize a vector b with random elements
    b = randi([1, 10], n, 1);

    % Calculate A1 and b1 for the system solution
    A1 = D - L * inv(D) * U - U * inv(D) * L - L * inv(D) * L - U * inv(D) * U;
    b1 = b + L * (inv(D) * b) + U * (inv(D) * b);

    % Solve the system
    x = A1 \ b1;

    % Compute A2 and A3 based on A1
    A2 = (diag(diag(A1)) + diag(diag(A1)) - A1) * inv(diag(diag(A1))) * A1;
    A3 = (diag(diag(A2)) + diag(diag(A2)) - A2) * inv(diag(diag(A2))) * A2;

    % Display results
    disp('Solution x:');
    disp(x);
    disp('Matrix A:');
    disp(A);
    disp('Matrix A1:');
    disp(A1);
    disp('Matrix A2:');
    disp(A2);
    disp('Matrix A3:');
    disp(A3);

    % Calculate and display norms
    norm_main_diag = norm(diag(D), 2);
    norm_sub_diag = norm(diag(L, -1), 2);
    norm_hyper_diag = norm(diag(U, 1), 2);
    fprintf('Norm of main diagonal: %f\n', norm_main_diag);
    fprintf('Norm of sub-diagonal: %f\n', norm_sub_diag);
    fprintf('Norm of hyper-diagonal: %f\n', norm_hyper_diag);

    % Simplified computational cost for Hadamard products
    % Assuming hypothetical direct Hadamard product operations
    cost_hadamard = n + 2 * (n - 1); % n for D.*b, (n-1) for L.*b and U.*b
    fprintf('Total computational cost estimate for Hadamard products: %d operations\n', cost_hadamard);
end
