function x = solve_triangular_system_experiment(A, B, epsilon)
    disp('Matrix A:');
    disp(A);
    disp(['Matrix Size: ', num2str(size(A))]);

    n = size(A, 1);
    x = zeros(n, size(B, 2));
    B_original = B; % Store the original B

    % Initialize arrays to store norms
    hyper_diag_norms = zeros(ceil(log2(n)), 1);
    sub_diag_norms = zeros(ceil(log2(n)), 1);

    for k = 1:log2(n)
        % Calculate norms of hyper- and sub-diagonals
        hyper_diag_norms(k) = norm(triu(A, 1), 'fro');
        sub_diag_norms(k) = norm(tril(A, -1), 'fro');

        % Regularized inverse of the diagonal
        D_inv = 1./(diag(A) + epsilon);

        % Update A and B
        A = D_inv .* (A + tril(A, -1) + triu(A, 1));
        B = D_inv .* (A + tril(A, -1) + triu(A, 1)) * B_original;

        % Update solution
        x = x + D_inv .* B;
    end

    % Display norms
    disp('Norms of hyper-diagonals:');
    disp(hyper_diag_norms);

    disp('Norms of sub-diagonals:');
    disp(sub_diag_norms);
    disp('---------------------------');
end