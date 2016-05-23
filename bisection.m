function [ x ] = bisection( f, a, b, n, epsilon )
%
% Inputs: f: An inline function.
c = f(a);
d = f(b);
if c*d > 0.0
    error('Please use other initial point.');
end

for i = 1:n
    x = (a + b)/2;
    y = f(x);
    if abs(y) < epsilon
        break;
    end
    fprintf('Iteration:%d: x = %f, y = %f\n', i, x, double(y));
    if c*y > 0.0
        a = x;
    else
        b = x;
    end
end

end

