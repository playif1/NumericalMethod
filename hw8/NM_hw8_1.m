%% Numerical Method Homework 8-1, cubic spline approximation, r04942044, Chien-Wen Huang.
% Set the input x and corresponding y, note that x and y should be column
% vectors with the same length.
x = (-10:0.25:10)';
xx = -10:0.25:10;
y = sin(x);

S = cubic_spline( x, y );
inner_points = 1;
plot_cubic_spline(x, S, inner_points);

% Compare with MATLAB spline.
yy = spline(x,y, xx);
plot(xx,yy, 'b-');

