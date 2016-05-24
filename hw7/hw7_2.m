%% Numerical Method Homework 7-2, R04942044, Jian-Wen Huang.
clc; close all; clear;

x_k = [1; 1; 1];
epsilon = 10e-7;
n = input('\nType the maximum times of iterations \n');
x_k(1) = input('\nType a initial x1 value\n');
x_k(2) = input('\nType a initial x2 value\n');
x_k(3) = input('\nType a initial x3 value\n');
x_k1 = x_k;

syms x1 x2 x3;
f1 = symfun(x1^2 + x2^2 + x3^2 - 1, [x1 x2 x3]);
f2 = symfun(x1^2 + x3^2 - 1/4, [x1 x2 x3]);
f3 = symfun(x1^2 + x2^2 - 4*x3, [x1 x2 x3]);
jac = jacobian([f1 f2 f3],[x1 x2 x3]);

for i = 1:n
    fprintf('Iteration %d, x1 = %f, x2 = %f, x3 = %f\n', i, double(x_k1(1)), double(x_k1(2)), double(x_k1(3)));    
    fx_k = [f1(x_k(1), x_k(2), x_k(3)); f2(x_k(1), x_k(2), x_k(3)); f3(x_k(1), x_k(2), x_k(3))];    
    x_k1 = x_k - jac(x_k(1), x_k(2), x_k(3))\fx_k;
    if double(norm(x_k1 - x_k)) < epsilon
        fprintf('Finished at iteration %d\n', i);
        break;
    end
    x_k = x_k1;
end
fprintf('x1 = %f, x2 = %f, x3 = %f\n', double(x_k1(1)), double(x_k1(2)), double(x_k1(3)));
