function plot_cubic_spline(x, S, inner_points)
% plots a cubic spline with break points x and coefficents S.s0, S.s1, S.s2, S.s3

n = length(x);

hold on;
for i=1:n-1
    xx = linspace(x(i),x(i+1),inner_points);
    xi = repmat(x(i),1,inner_points);
    yy = S.s0(i) + S.s1(i)*(xx-xi) + S.s2(i)*(xx-xi).^2 + S.s3(i)*(xx - xi).^3;
    plot(xx,yy,'ro')
    plot(x(i),0,'r*');
end

