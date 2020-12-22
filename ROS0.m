function u = ROS0(h, x_n, tau, t_n, mu, mu_1, mu_2, k, f)

%создание сеток
x = 0:h:x_n;
t = 0:tau:t_n;
n_x = length(x);
n_t = length(t);

%создание матрицы лямбда
Lambda = zeros(n_x);
for i = 2:n_x-1
    Lambda(i,i) = -2*k/(h^2);
    Lambda(i,i-1) = k/(h^2);
    Lambda(i,i+1) = k/(h^2);
end

%краевые условия
u(:,1) = mu_1(t);
u(:,n_x) = mu_2(t);

%начальное условие
u(1,1:n_x) = mu(x(1:n_x));

L = Lambda(2:n_x-1,:);
for i = 2:n_t
    f_temp = double(f(x(2:n_x-1),t(i)));
    u_temp = (L*(u(i-1,:))')';
    u(i,2:n_x-1) =  u(i-1,2:n_x-1)+ tau*u_temp + tau*f_temp;
end
end