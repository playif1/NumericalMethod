function [ x ] = newton( f, x0, n, epsilon, xStar )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
def = diff(f);
x = x0;
cr = 1;

for i = 1:n
    x0 = x;
    x = x0 - f(x0)/def(x0);
    y = f(x);
    if i > 1
        cr = abs(x-xStar)/abs(x0-xStar);
    end    
    fprintf('Iteration:%d: x = %f, y = %f\n', i, double(x), double(y));
    fprintf('The convergence rate = %f\n', double(cr));
    if abs(y) < epsilon
        break;
    end
end

end

