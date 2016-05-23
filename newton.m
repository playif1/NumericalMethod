function [ x ] = newton( f, x0, n, epsilon )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
def = diff(f);
x = x0;
for i = 1:n
    x0 = x;
    x = x0 - f(x0)/def(x0);
    y = f(x);
    fprintf('Iteration:%d: x = %f, y = %f\n', i, double(x), double(y));
    if abs(y) < epsilon
        break;
    end
end

end

