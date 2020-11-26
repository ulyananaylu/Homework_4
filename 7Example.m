clear; clc; close('all');

%% Inputs:
% specify step
h = 0.01;

% specify boundaries
a = -1;
b = 1;

% specify amount of coordinate functions
n = 3;

%% Calculations:
% Using Galerkin method

y = Galerkin(n, a, b);

disp(y);

X = (a:h:b)';

% Y = zeros(length(X),1);
% 
% for i = 1:length(X)
%     x = X(i);
%     Y(i) = subs(y);
% end


figure;
fplot(y, [a b], 'b');
title('Result plot');
hold('on');

X = linspace(a,b,length(X));
solinit = bvpinit(X, [0 0]);
sol = bvp4c(@fun, @bc2, solinit);
Yex = deval(sol, X);
plot(X, Yex(1,:), 'r');
legend('Me', 'Matlab');

dy = diff(y);
ddy = diff(dy);

dYex = diff([X Yex(1,:)]);
ddYex = diff([X Yex(1,:)], 2);

fprintf('%10s%30s%30s%30s\n', 'Xi', 'my', 'matlab');
fprintf('\n');
for i = 1:length(X)
    mydL = P*ddy + Q*dy + R*y - F;
    x = X(i);
    dL = subs(mydL);
    
    matlabdL = p(X(i))*ddYex(i) + q(X(i))*dYex(i) + r(X(i))*Yex(i) - f(X(i));
    
    fprintf('%10.7f%30.7f%30.7f%30.7f\n', X(i), dL, matlabdL);
    fprintf('\n');
end

% fprintf('%10s%30s%30s%30s\n', 'Xi', 'Yi', 'Yexi', 'abs(Yi - Yexi)');
% for i = 1:length(x)
%     fprintf('%10.7f%30.7f%30.7f%30.7f\n', X(i), Y(i), Yex(1,i), abs(Y(i) - Yex(1,i)));
% end

function tmp = p(x)
    tmp = (x - 2) / (x + 2);
end

function tmp = q(x)
    tmp = x;
end

function tmp = r(x)
    tmp = 1 - sin(x);
end

function tmp = f(x)
    tmp = x.^2;
end

function F=fun(x,y)
F=[
    y(2);
%     (3 -x) * (2 - x - exp(x/2) * y(1) - (1 + x/2) * y(2))
    (x*x -x * y(2) - (1 - sin(x))* y(1)) *(x+2)/(x-2)
   ];
end

function res=bc2(ya, yb)
    res = [
            ya(1);
            yb(1)
           ];
end
