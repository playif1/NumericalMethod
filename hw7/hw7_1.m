%% Numerical Method Homework 7-1, R04942044, Jian-Wen Huang.
clc; close all; clear;

% Input the arguments for the bisection method and the Newton's method.
syms x;
Fun = input('Type a function \n');
a = input('\nType a initial left point for the Bisection method \n');
b = input('\nType a initial right point for the Bisection method \n');
n = input('\nType the maximum times of iterations \n');
epsilon = input('\nType the epsilon value for quicker convergence \n');
x0 = input('\nType a initial x value for the Newton''s method\n');

f = symfun(Fun, x);
plot(a:0.1:b, f(a:0.1:b));
xStar = solve(f);
xStar = xStar(xStar<b & xStar>a);
% Bisection Method
fprintf('\nThe Bisection Method:\n');
xb = bisection(f, a, b, n, epsilon, xStar);
fprintf('\nAnswer by Bisection method = %f\n', double(xb));

% Newton's Method
fprintf('\nThe Newton''s Method:\n');
xn = newton(f, x0, n, epsilon, xStar);
fprintf('\nAnswer by Newton''s method = %f\n', double(xn));