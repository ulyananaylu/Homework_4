clear; clc; close('all');

%% Inputs:
% specify p, q, r, f p*y'' + q*y' + r*y = f in respective functions
% One note: in this method assumed that L = p*y'' + q*y' + r*y = f. Thanks Vlad for catching this)

% specify step
h = 0.001;

% specify boundaries
a = -1;
b = 1;

% specify boundary conditions
alpha1 = 1;
alpha2 = 0;
alpha = 0;

beta1 = 1;
beta2 = 0;
beta = 0;

[X, Y] = Tridiagonal(a, b, h, alpha1, alpha2, alpha, beta1, beta2, beta);

figure;
plot(X, Y, 'b');
title('Result plot');
hold('on');

x = linspace(a,b,length(X));
solinit = bvpinit(x, [0 0]);
sol = bvp4c(@fun, @bc2, solinit);
Yex = deval(sol, x);
plot(x, Yex(1,:), 'r');
legend('Me', 'Matlab');

fprintf('%10s%30s%30s%30s\n', 'Xi', 'Yi', 'Yexi', 'abs(Yi - Yexi)');
for i = 1:length(x)
    fprintf('%10.7f%30.7f%30.7f%30.7f\n', X(i), Y(i), Yex(1,i), abs(Y(i) - Yex(1,i)));
end


function F=fun(x,y)
F=[
    y(2);
%     (3-x) * (2 - x - exp(x/2) * y(1) - (1 + x/2) * y(2))
    (x*x -x * y(2) - (1 - sin(x))* y(1)) *(x+2)/(x-2)
   ];
end

function res=bc2(ya, yb)
    res = [
            ya(1);
            yb(1)
           ];
end
