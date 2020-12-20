global fun; 
clc
format long;

k = 1;

tre = solex();
check = tre(1, 2:4);

    n = 7;
    
    [func, dfunc, ddfunc] = coord_func_and_deriv(k, n);
    points = getChebZeros(n, -1, 1);
    
    LW_oper = sym(zeros(1, n));
    for i = 1:n
        LW_oper(i) = my_operator(func(i), dfunc(i), ddfunc(i));
    end
    
    LWW_at_point = values_operator(LW_oper, points, n);
    
    syms x;
    f = 1 + x;
    FW_at_point = zeros(1, n);
    fun = f;
    FW_at_point = funcX(points);
    
    c_j = linsolve(LWW_at_point, FW_at_point');
    disp(n);
    c_j
    FW_at_point

    
    syms y_n;
    y_n = sum(c_j .* func);
    
    fun = y_n;
    
    X = -1:(2/100):1;
    
    Y = funcX(X);
    
    plot(X,Y, 'green');
    hold on;
   
function tmp = myjacobi(k, n)
    pj = sym(zeros(n + 1, 1));
    syms x;
    pj(1) = 1;
    switch n
        case 0
            tmp = pj;
        case 1
            pj(2) = (k + 1) * x;
            tmp = pj;
        otherwise
            pj(2) = (k + 1) * x;
            for i = 2:n
                pj(i + 1) = ((i + k)/(i + 2 * k)) * (2 + (2 * k - 1) / i) * x * pj(i) - ...
                    ((i + k) / (i + 2 * k)) * (1 + (k - 1) / i) * pj(i - 1);
            end
            tmp = pj;
    end
end

function [phi, dphi, ddphi] = coord_func_and_deriv(k, n)
    phi = sym(zeros(n, 1));
    dphi = sym(zeros(n, 1));
    ddphi = sym(zeros(n, 1));
    syms x;
    jac = myjacobi(k, n - 1);
    djac = myjacobi(k - 1, n);
    for i = 1:n
        phi(i) = (1 - x^2) * jac(i);
        phi(i) = collect(phi(i));
        dphi(i) = -2 * (i) * (1 - x^2)^(k - 1) * djac(i + 1);
        dphi(i) = collect(dphi(i));
        ddphi(i) = -2 * (i) * ((k - 1) * (1 - x^2)^(k - 2) * (-2 * x) * djac(i + 1) + ...
            (1 - x^2)^(k - 1) * (i + 2*(k - 1) + 1) / 2) * jac(i);
        ddphi(i)=collect(ddphi(i));
    end
end

function oper = my_operator(coord, dcoord, ddcoord)
syms oper;
syms koef;
syms x;
koef = [-1/(2 + x) , 1  / (2+x)^2, cos(x)];
oper = koef(1) * ddcoord + koef(2) * dcoord + koef(3) * coord;
end

function mat = values_operator(operator, points, n)
global fun;
syms temp;
mat = zeros(n);

for i = 1:n
    for j = 1:n
        fun = operator(j);
        mat(i, j) = funcX(points(i));
    end
end
end

function y = funcX(dot)
global fun;
syms x;
y = subs(fun, x, dot);
y = double(y);
end

function t = getChebZeros(n,a,b)
t=zeros(1,n);
for i=1:n
    t(i)=(a+b)/2+(b-a)/2*cos((2*i-1)/n/2*pi);
end
end






