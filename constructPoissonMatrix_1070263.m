function A = constructPoissonMatrix_1070263(X, ~, h)
    % Construct the Poisson matrix A for the given grid and boundary conditions
    n = size(X, 1);
    m = size(X, 2);
    A = sparse(n * m, n * m);

    for i = 1:n
        for j = 1:m
            idx = (i - 1) * m + j;

            % Diagonal element
            A(idx, idx) = -4 / h^2;

            % Neighbors
            if i > 1
                A(idx, idx - m) = 1 / h^2;
            end
            if i < n
                A(idx, idx + m) = 1 / h^2;
            end
            if j > 1
                A(idx, idx - 1) = 1 / h^2;
            end
            if j < m
                A(idx, idx + 1) = 1 / h^2;
            end
        end
    end

    % Apply boundary conditions
    A(1:m, :) = 0;
    A(end - m + 1:end, :) = 0;
    A(:, 1:m) = 0;
    A(:, end - m + 1:end) = 0;
end


