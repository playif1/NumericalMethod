%% Numerical Method Homework 7-1, R04942044, Jian-Wen Huang.
syms x;
a = 1.5; 
b = 4;
n = 200;
epsilon = 10e-6;

Fun = input('\nType a function \n');
x0 = input('\nType a initial x value\n');

f = symfun(Fun, x);

% Bisection Method
xb = bisection(f, a, b, n, epsilon);
fprintf('Answer by Bisection method = %f\n', double(xb));
plot(a:0.1:b, f(a:0.1:b));

% Newton's Method
xn = newton(f, x0, n, epsilon);
fprintf('Answer by Newton''s method = %f\n', double(xn));