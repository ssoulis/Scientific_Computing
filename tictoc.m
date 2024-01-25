% Number of runs for benchmarking
num_runs = 10;

% Benchmarking without preconditioning
time_noprec = zeros(num_runs, 1);
for run = 1:num_runs
    tic;
    [x_noprec, flag_noprec, relres_noprec, iter_noprec, resvec_noprec, errvec_noprec] = pcg_myid_1070263(A, b, 1e-6, maxiter, 'none', [], xsol);
    time_noprec(run) = toc;
end

% Benchmarking with incomplete Cholesky preconditioning
time_ic = zeros(num_runs, 1);
for run = 1:num_runs
    tic;
    [x_ic, flag_ic, relres_ic, iter_ic, resvec_ic, errvec_ic] = pcg_myid_1070263(A, b, 1e-6, maxiter, 'ichol', [], xsol);
    time_ic(run) = toc;
end

% Benchmarking with custom preconditioning
time_custom = zeros(num_runs, 1);
for run = 1:num_runs
    tic;
    [x_custom, flag_custom, relres_custom, iter_custom, resvec_custom, errvec_custom] = pcg_myid_1070263(A, b, 1e-6, maxiter, 'custom', [], xsol);
    time_custom(run) = toc;
end

% Display average execution times
avg_time_noprec = mean(time_noprec);
avg_time_ic = mean(time_ic);
avg_time_custom = mean(time_custom);

fprintf('Average Execution Time (No Preconditioning): %.6f seconds\n', avg_time_noprec);
fprintf('Average Execution Time (IC(0) Preconditioning): %.6f seconds\n', avg_time_ic);
fprintf('Average Execution Time (Custom Preconditioning): %.6f seconds\n', avg_time_custom);
% Initialize tables
table_noprec = table();
table_ic = table();
table_custom = table();

% Without Preconditioning
table_noprec.Method = {'No Preconditioning'};
table_noprec.Flag = flag_noprec;
table_noprec.Final_RelResidual = relres_noprec(end);
table_noprec.Final_RelError = errvec_noprec(end);
table_noprec.Avg_ExecutionTime = mean(exec_times_noprec);  % Average execution time

% With Incomplete Cholesky Preconditioning
table_ic.Method = {'IC(0) Preconditioning'};
table_ic.Flag = flag_ic;
table_ic.Final_RelResidual = relres_ic(end);
table_ic.Final_RelError = errvec_ic(end);
table_ic.Avg_ExecutionTime = mean(exec_times_ic);  % Average execution time

% With Custom Preconditioning
table_custom.Method = {'Custom Preconditioning'};
table_custom.Flag = flag_custom;
table_custom.Final_RelResidual = relres_custom(end);
table_custom.Final_RelError = errvec_custom(end);
table_custom.Avg_ExecutionTime = mean(exec_times_custom);  % Average execution time

% Display tables
disp('Table Without Preconditioning:');
disp(table_noprec);

disp('Table With Incomplete Cholesky Preconditioning:');
disp(table_ic);

disp('Table With Custom Preconditioning:');
disp(table_custom);