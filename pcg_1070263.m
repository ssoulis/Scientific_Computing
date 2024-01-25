function [x, flag, relres, iter, resvec, errvec] = pcg_1070263(A, b, tol, maxit, preconditioner, x0, xsol, varargin)
    % PCG_MYID_1070263   Modified Preconditioned Conjugate Gradients Method.
    %   X = PCG_MYID_1070263(A, B, TOL, MAXIT, PRECONDITIONER, X0, XSOL) attempts to solve the
    %   system of linear equations A*X=B for X with the specified preconditioner.
    %   It returns the A-norm of the error vector.
    %
    %   Other input and output arguments are the same as the original PCG
    %   function.

    if nargin < 2
        error(message('MATLAB:pcg:NotEnoughInputs'));
    end

    % Determine whether A is a matrix or a function.
    if isnumeric(A) || islogical(A)
        atype = 'matrix';
        afun = [];
    else
        atype = 'function';
        afun = A;
    end

    if strcmp(atype, 'matrix')
        % Check matrix and right-hand side vector inputs have appropriate sizes
        [m, n] = size(A);
        if (m ~= n)
            error(message('MATLAB:pcg:NonSquareMatrix'));
        end
        if ~isequal(size(b), [m, 1])
            error(message('MATLAB:pcg:RSHsizeMatchCoeffMatrix', m));
        end
    else
        m = size(b, 1);
        n = m;
        if ~iscolumn(b)
            error(message('MATLAB:pcg:RSHnotColumn'));
        end
    end

    % Assign default values to unspecified parameters
    if (nargin < 3) || isempty(tol)
        tol = 1e-6;
    end
    warned = 0;
    if tol <= eps
        warning(message('MATLAB:pcg:tooSmallTolerance'));
        warned = 1;
        tol = eps;
    elseif tol >= 1
        warning(message('MATLAB:pcg:tooBigTolerance'));
        warned = 1;
        tol = 1 - eps;
    end
    if (nargin < 4) || isempty(maxit)
        maxit = min(n, 20);
    end
    maxit = max(maxit, 0);

    % Check for all zero right hand side vector => all zero solution
    n2b = norm(b);                     % Norm of rhs vector, b
    if (n2b == 0)                      % if    rhs vector is all zeros
        x = zeros(n, 1);                % then  solution is all zeros
        flag = 0;                      % a valid solution has been obtained
        relres = 0;                    % the relative residual is actually 0/0
        iter = 0;                      % no iterations need be performed
        resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
        errvec = 0;                    % Initialize the error vector
        if (nargout < 2)
            itermsg('pcg', tol, maxit, 0, flag, iter, NaN);
        end
        return
    end

    if ((nargin >= 5) && ~isempty(preconditioner))
        existM1 = 1;
        [m1type, m1fun] = iterchk(preconditioner);
        if strcmp(m1type, 'matrix')
            if ~isequal(size(preconditioner), [m, m])
                error(message('MATLAB:pcg:WrongPrecondSize', m));
            end
        end
    else
        existM1 = 0;
        m1type = 'matrix';
    end

    if ((nargin >= 6) && ~isempty(x0))
        if ~isequal(size(x0), [n, 1])
            error(message('MATLAB:pcg:WrongInitGuessSize', n));
        else
            x = x0;
        end
    else
        x = zeros(n, 1);
    end

       % Set up for the method
    flag = 1;
    xmin = x;                          % Iterate which has minimal residual so far
    imin = 0;                          % Iteration at which xmin was computed
    tolb = tol * n2b;                  % Relative tolerance
    r = b - iterapp('mtimes', A, 'matrix', x, varargin{:});
    normr = norm(r);                   % Norm of residual
    normr_act = normr;

    % Apply preconditioner
    switch preconditioner
        case 'ichol'
            [L, U] = incompleteCholesky(A);
        case 'custom'
            [L, U] = customPreconditioner(A);
        otherwise
            L = speye(size(A));
            U = speye(size(A));
    end

    r = U \ (L \ r);
    
    if (normr <= tolb)                 % Initial guess is a good enough solution
        flag = 0;
        relres = normr / n2b;
        iter = 0;
        resvec = normr;
        errvec = norm(xsol - x) / norm(xsol);
        if (nargout < 2)
            itermsg('pcg', tol, maxit, 0, flag, iter, relres);
        end
        return
    end

    resvec = zeros(maxit + 1, 1);         % Preallocate vector for norm of residuals
    resvec(1) = normr;               % resvec(1) = norm(b-A*x0)
    normrmin = normr;                  % Norm of minimum residual
    rho = 1;
    stag = 0;                          % stagnation of the method
    moresteps = 0;
    maxmsteps = min([floor(n/50), 5, n - maxit]);
    maxstagsteps = 3;

    % loop over maxit iterations (unless convergence or failure)
    for ii = 1:maxit
        if existM1
            y = iterapp('mldivide', m1fun, m1type, r, varargin{:});
            if ~allfinite(y)
                flag = 2;
                break
            end
        else % no preconditioner
            y = r;
        end

        rho1 = rho;
        rho = r' * y;
        if ((rho == 0) || isinf(rho))
            flag = 4;
            break
        end
        if (ii == 1)
            p = y;
        else
            beta = rho / rho1;
            if ((beta == 0) || isinf(beta))
                flag = 4;
                break
            end
            p = y + beta * p;
        end
        q = iterapp('mtimes', afun, atype, p, varargin{:});
        pq = p' * q;
        if ((pq <= 0) || isinf(pq))
            flag = 4;
            break
        else
            alpha = rho / pq;
        end
        if isinf(alpha)
            flag = 4;
            break
        end

        % check for convergence
        if (normr_act <= tolb || stag >= maxstagsteps || moresteps)
            r = b - iterapp('mtimes',afun,atype,x,varargin{:});
            normr_act = norm(r);
            resvec(ii+1, 1) = normr_act;
            errvec(ii+1, 1) = norm(xsol - x, A) / norm(xsol, A);  % A-norm of the error

            if (normr_act <= tolb)
                flag = 0;
                iter = ii;
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('MATLAB:pcg:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = ii;
                    break;
                end
            end
        end

        % Check for stagnation of the method    
        if (norm(p) * abs(alpha) < eps * norm(x))
            stag = stag + 1;
        else
            stag = 0;
        end

        x = x + alpha * p;             % form new iterate
        r = r - alpha * q;
        normr = norm(r);
        normr_act = normr;
        resvec(ii + 1) = normr;

        % check for convergence
        if (normr_act <= tolb || stag >= maxstagsteps || moresteps)
            r = b - iterapp('mtimes', afun, atype, x, varargin{:});
            normr_act = norm(r);
            resvec(ii + 1) = normr_act;
            if (normr_act <= tolb)
                flag = 0;
                iter = ii;
                errvec = norm(xsol - x) / norm(xsol);
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('MATLAB:pcg:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = ii;
                    errvec = norm(xsol - x) / norm(xsol);
                    break;
                end
            end
        end

        if (normr_act < normrmin)      % update minimal norm quantities
            normrmin = normr_act;
            xmin = x;
            imin = ii;
        end
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
    end  % for ii = 1:maxit

    if isempty(ii)
        ii = 0;
    end

    % returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else
        r_comp = b - iterapp('mtimes', afun, atype, xmin, varargin{:});
        if norm(r_comp) <= normr_act
            x = xmin;
            iter = imin;
            relres = norm(r_comp) / n2b;
            errvec = norm(xsol - x) / norm(xsol);
        else
            iter = ii;
            relres = normr_act / n2b;
            errvec = norm(xsol - x) / norm(xsol);
        end
    end

    % truncate the zeros from resvec
    if ((flag <= 1) || (flag == 3))
        resvec = resvec(1:ii + 1);
        errvec = errvec(1:ii + 1);
    else
        resvec = resvec(1:ii);
        errvec = errvec(1:ii);
    end

    % only display a message if the output flag is not used
    if (nargout < 2)
        itermsg('pcg', tol, maxit, ii, flag, iter, relres);
    end

end  % pcg_myid_1070263
