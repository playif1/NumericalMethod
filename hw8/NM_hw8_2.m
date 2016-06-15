%% Numerical Method Homework 8-2, quadratic least square approximation, r04942044, Chien-Wen Huang.
x = (-10:0.25:10)';
y = sin(x);
[a, b, c] = quadratic_least_square(x, y);
plot(x, a*x.^2 + b*x + c, 'r-', x, y, 'b-');

fprintf('the total error for quatratic least square approximation is %f\n', sum((y - (a*x.^2 + b*x + c)).^2));