function [ x ] = bisection( f, a, b, n, epsilon, xStar )
%
% Inputs: f: An inline function.
c = f(a);
d = f(b);
if c*d > 0.0
    error('Please use other initial points.');
end
x = (a + b)/2;
cr = 1;

for i = 1:n
    x0 = x;
    x = (a + b)/2;
    y = f(x);
    if abs(y) < epsilon
        break;
    end
    if i > 1
        cr = abs(x-xStar)/abs(x0-xStar);
    end
    fprintf('Iteration:%d: x = %f, y = %f\n', i, x, double(y));
    fprintf('The convergence rate = %f\n', double(cr));
    if c*y > 0.0
        a = x;
    else
        b = x;
    end
end

end

