m = 5;
x= [0:1:2*m-1]'/m; 
y = x.^4 + x.^3 - 2*x.^2 + 4*x + log(x+1) - 3*cos(x);
z = pi*(x-1);
n=5;
a = []; b = [];
for k=0:n
a = [a ; 1/m*sum(y.*cos(k*z))];
end
for k=1:n-1
b = [b; 1/m*sum(y.*sin(k*z))] ;
end
approx=a(1)/2 + cos(z*(1:n))*a(2:n+1) + sin(z*(1:n-1))*b(1:n-1);
% z*(1:n) is a matrix
[y approx]
E = sum((y-approx).^2)
plot(x,y,'b-',x,approx,'r--');