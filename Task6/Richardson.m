function [u, x, d] = Richardson(x_0, x_n, q, r, f, alpha, beta,epsilon)

%x_0 - первый узел, x_n - последний узел, q r f - из дифура
%alpha, beta - отвечают за краевые условия, epsilon - заданная точность
    
    %начальный шаг
    h = (x_n - x_0)/10;
    
    %применение метода для получения значения на первой сетке
    [u_1, x_1] = net_met(x_0, h, x_n, q, r, f, alpha, beta);
    
    d = 1;%погрешность
    iter = 1;%количество итераций
    
    %условие точности
    while max(abs(d))>epsilon
        h = h/2;%сгущаем сетку
        [u_2, x_2] = net_met(x_0, h, x_n, q, r, f, alpha, beta);%вычисляем значения на новой сетке
        iter = iter +1;
        dd = zeros(3,1);%матрица с тремя погрешностями
        
        %нахождение максимальной погрешности между сетками
        dd(1) = (u_2(1) - u_1(1))/(3);
        m = 0;
        for j = 2:length(x_1)
            dd(3) = (u_2(2*j-1) -u_1(j))/(3);
            dd(2) = (dd(1) + dd(3))/2;
            for i = 1:3
                if max(abs(dd(i)))>m
                    m = max(abs(dd(i)));
                    d = dd(i);
                end
            end
            dd(1) = dd(3);
        end
        
        u_1 = u_2; x_1 = x_2;
    end
    disp('Количество итераций:');
    disp(iter);
    u = u_2;
    x = x_2;
end
