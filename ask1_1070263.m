clear; clc;

% Step 1: Generate data
n = 100 * 2.^(0:6);
for j = 1:numel(n)
    A = randn(n(j));
    B = randn(n(j));
    u = randn(n(j), 1);
    v = randn(n(j), 1);
    fn_ger = @() A + B * (u * v');
    t(j) = timeit(fn_ger);
end

% Step 2: Use polyfit with weights to fit a cubic function
degree = 3;
[coefficients, S, mu] = polyfit(n, t, degree);

% Step 3: Extract coefficients
alpha3 = coefficients(1);
alpha2 = coefficients(2);
alpha1 = coefficients(3);
alpha0 = coefficients(4);

% Display the coefficients
disp('Cubic function coefficients:');
disp(['alpha3: ', num2str(alpha3)]);
disp(['alpha2: ', num2str(alpha2)]);
disp(['alpha1: ', num2str(alpha1)]);
disp(['alpha0: ', num2str(alpha0)]);

% Plot the data and the fitted cubic function
figure;
plot(n, t, '-x', 'DisplayName', 'Measured Data');
hold on;
fitted_curve = polyval(coefficients, n, S, mu);
plot(n, fitted_curve, 'r-', 'DisplayName', 'Fitted Cubic Function');
xlabel('n');
ylabel('Time (sec)');
title('Fitting Cubic Function to Execution Times');
legend('show');

% Given coefficients from polyfit
alpha3 = coefficients(1);
alpha2 = coefficients(2);
alpha1 = coefficients(3);
alpha0 = coefficients(4);

% Calculate predicted times for n in [100:100:2000]
n_range1 = 100:100:2000;
predicted_times_range1 = polyval(coefficients, n_range1, S, mu);

% Calculate predicted times for n in [150:100:1550]
n_range2 = 150:100:1550;
predicted_times_range2 = polyval(coefficients, n_range2, S, mu);

% Display the predicted times for both ranges
disp('Predicted execution times for n in [100:100:2000]:');
disp(predicted_times_range1);

disp('Predicted execution times for n in [150:100:1550]:');
disp(predicted_times_range2);

% Given coefficients from polyfit
alpha3 = coefficients(1);
alpha2 = coefficients(2);
alpha1 = coefficients(3);
alpha0 = coefficients(4);

% Calculate predicted times for n in [100:100:2000]
n_range1 = 100:100:2000;
predicted_times_range1 = polyval(coefficients, n_range1, S, mu);

% Calculate predicted times for n in [150:100:1550]
n_range2 = 150:100:1550;
predicted_times_range2 = polyval(coefficients, n_range2, S, mu);

% Measure actual execution times for n in [100:100:2000]
actual_times_range1 = zeros(size(n_range1));
for j = 1:numel(n_range1)
    A = randn(n_range1(j));
    B = randn(n_range1(j));
    u = randn(n_range1(j), 1);
    v = randn(n_range1(j), 1);
    fn_ger = @() A + B * (u * v');
    actual_times_range1(j) = timeit(fn_ger);
end

% Measure actual execution times for n in [150:100:1550]
actual_times_range2 = zeros(size(n_range2));
for j = 1:numel(n_range2)
    A = randn(n_range2(j));
    B = randn(n_range2(j));
    u = randn(n_range2(j), 1);
    v = randn(n_range2(j), 1);
    fn_ger = @() A + B * (u * v');
    actual_times_range2(j) = timeit(fn_ger);
end

% Plotting results
figure;

% Plot actual execution times for n in [100:100:2000]
subplot(2, 1, 1);
plot(n_range1, actual_times_range1, 'o-', 'DisplayName', 'Actual Times');
hold on;
plot(n_range1, predicted_times_range1, 'r-', 'DisplayName', 'Predicted Times');
xlabel('n');
ylabel('Time (sec)');
title('Accuracy Test - n in [100:100:2000]');
legend('show');

% Plot actual execution times for n in [150:100:1550]
subplot(2, 1, 2);
plot(n_range2, actual_times_range2, 'o-', 'DisplayName', 'Actual Times');
hold on;
plot(n_range2, predicted_times_range2, 'r-', 'DisplayName', 'Predicted Times');
xlabel('n');
ylabel('Time (sec)');
title('Accuracy Test - n in [150:100:1550]');
legend('show');

% Adjust layout
subplot(2, 1, 1);
ylim([0 max([actual_times_range1, predicted_times_range1])]);

% Adjust layout
subplot(2, 1, 2);
ylim([0 max([actual_times_range2, predicted_times_range2])]);

