function [ a, b, c ] = quadratic_least_square( x, y )
% x, y should be column vectors.
%   Detailed explanation goes here
Sx = sum(x);
Sxx = sum(x.^2);
Sxxx = sum(x.^3);
Sxxxx = sum(x.^4);
Sy = sum(y);
Sxy = sum(x.*y);
Sxxy = sum(y.*(x.^2));

a = ((Sxy-Sx*Sy)*(Sxxx-Sxx*Sxx) - (Sxxy-Sxx*Sy)*(Sxx-Sx*Sx)) / ((Sxxx-Sxx*Sx)*(Sxxx-Sxx*Sxx) - (Sxxxx-Sxx*Sxx)*(Sxx-Sx*Sx));
b = ((Sxxy-Sy*Sxx)*(Sxxx-Sxx*Sx) - (Sxy-Sx*Sy)*(Sxxxx-Sxx*Sxx)) / ((Sxx-Sx*Sx)*(Sxxxx-Sxx*Sxx) - (Sxxx-Sxx*Sxx)*(Sxxx-Sxx*Sx));
c = ((Sxxy*Sxx-Sy*Sxxxx)*(Sxxx*Sxxx-Sxx*Sxxxx) - (Sxxy*Sxxx-Sxy*Sxxxx)*(Sxxx*Sxx-Sx*Sxxxx)) / ((Sxx*Sxxx-Sx*Sxxxx)*(Sxxx*Sxx-Sx*Sxxxx) - (Sxx*Sxx-Sxxxx)*(Sxxx*Sxxx-Sxx*Sxxxx));
end

