clear; clc;

% Generate data for n = [100:100:2000]
n = 100:100:2000;
for j = 1:numel(n)
    A = randn(n(j));
    B = randn(n(j));
    u = randn(n(j), 1);
    v = randn(n(j), 1);
    fn_ger = @() A + B * (u * v');
    t(j) = timeit(fn_ger);
end

% Fit polynomials of degree 2, 3 (cubic), and 4 with weights
degree2_coefficients = polyfit(n, t, 2);
degree4_coefficients = polyfit(n, t, 4);

% Fit cubic polynomial with weights
[coefficients, S, mu] = polyfit(n, t, 3);

% Calculate predicted times for all n values
predicted_times_cubic = polyval(coefficients, n, S, mu);
predicted_times_degree2 = polyval(degree2_coefficients, n);
predicted_times_degree4 = polyval(degree4_coefficients, n);

% Plotting results
figure;

% Plot actual execution times and predicted times for all n
plot(n, t, 'o-', 'DisplayName', 'Actual Times');
hold on;
plot(n, predicted_times_cubic, 'r-', 'DisplayName', 'Cubic Prediction');
plot(n, predicted_times_degree2, 'g-', 'DisplayName', 'Degree 2 Prediction');
plot(n, predicted_times_degree4, 'b-', 'DisplayName', 'Degree 4 Prediction');
xlabel('n');
ylabel('Time (sec)');
title('Comparison of Polynomial Predictions');
legend('show');
