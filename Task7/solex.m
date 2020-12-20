function sol = solex()
%Задано уравнение вида p_eq(x)y''+q_eq(x)*y'+r_eq(x)*y=f_eq(x)
%Требуется привести к виду y''+q(x)*y'+r(x)*y=f(x)
%Объявляем коэффициенты исходного уравнения



    function [X] = p_eq(x)
         X = -1 / (2 + x);
    end
    function [X] = q_eq(x)
        X = 1 / (2 + x)^2;
    end
    function [X] = r_eq(x)
        X = cos(x);
    end
    function [X] = f_eq(x)
        X = 1 + x;
    end
%Строим коэффициенты нового уравнения
    function [X] = q(x)
        X = q_eq(x) / p_eq(x);
        
    end
    function [X] = r(x)
        X= r_eq(x) / p_eq(x);
    end
    function [X] = f(x)
        X = f_eq(x) / p_eq(x);
    end

clc;
format long;
%Задаем промежуток, на котором решается уравнение
a = -1;
b = 1;
%Задаем коэффициенты в граничных условиях
%   alpha1*u(a)-alpha2*u'(a)=alpha
%   beta1*u(b)+beta2*u'(b)=beta
alpha1 = 1;
alpha2 = 0;
alpha = 0;
beta1 = 1;
beta2 = 0;
beta = 0;  

%solution solver bvp4c()
%   y1'=y2;
%   y2'=f(x)-q(x)*y2-r(x)*y1

%Задаем правые части уравнений
    function tmp = twoode(x,y)
        tmp = [ y(2)
            f(x)-q(x)*y(2)-r(x)*y(1)];
    end
%Граничные условия следует выписать так, чтобы в правых частях стояли  нули
    function res = twobc(ya,yb)
        res = [ alpha1*ya(1)-alpha2*ya(2)-alpha
            beta1*yb(1)+beta2*yb(2)-beta];
    end
n = 100;
solinit = bvpinit(linspace(a, b, n+1),[0.6 0]);%Выбор начальной сетки и начального приближения к решению
options = bvpset('Reltol', 1e-06, 'abstol', 1e-06, 'Nmax', 1e10); %здесь задается точность решения - оценивается невязка
sol = bvp4c(@twoode, @twobc, solinit, options);%Собственно обращение к солверу

x = linspace(a, b, n+1); %Задаем сетку,в узлах которой требуется получить решение
ExSol = deval(sol, x);  %Интерполяция - вычисление решения в заданных точках сетки
sol = ExSol;
%disp('ExSol');
%disp(ExSol(1,:));
mas = [x; ExSol(1,:)];
fprintf('%5.2f  %10.5f \n', mas);
hold on;
plot(x, ExSol(1,:),'red');
hold on;

end
