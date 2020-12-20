global fun; 
clc
format long;

k = 1;

tre = solex(); 
check = tre(1, 2:4);


n = 7;

[func, dfunc, ddfunc] = coord_func_and_deriv(k, n);
    
LW_oper = sym(zeros(1, n)); 

    LWW_scalar = scalar_operator(func, dfunc, n); 
    
  
    syms x;
    f = 1 + x;
    FW_scalar = zeros(1, n);
    for i = 1:n
        fun = f * func(i);
        FW_scalar(i) = integral(@funcX, -1, 1);
    end
    
    c_j = linsolve(LWW_scalar, FW_scalar');
    disp(n);
    c_j
    LWW_scalar
    FW_scalar   

    
    syms y_n;
    y_n = sum(c_j .* func);
    
    fun = y_n;
    
    X = -1:(2/100):1;
    
    Y = funcX(X);
    
    plot(X,Y, 'blue');
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


function mat = scalar_operator(func, dfunc, n)
global fun;
syms temp;
mat = zeros(n);

syms x;
p = 1 / (2 + x);
r = cos(x);

for i = 1:n
    for j = 1:n
        
        
        
        temp =  p *  dfunc(i) * dfunc(j) + r * func(i) * func(j);

        fun = temp;
        mat(i, j) = integral(@funcX, -1, 1);
    end
end
end

function y = funcX(dot)
global fun;
syms x;
y = subs(fun, x, dot);
y = double(y);
end










