function u = CROS(h, x_n, tau, t_n, mu, mu_1, mu_2, k, f)

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
u(1,:) = mu(x);

%L = Lambda(2:n_x-1,:);
E = eye(size(Lambda));
Lambda(1,1) = 1;
Lambda(n_x,n_x) = 1;
for i = 2:n_t
    b = zeros(n_x,1);
    b(1) = u(i,1);
    b(n_x) = u(i,n_x);
    f_temp = double(f(x,t(i)));
    A = (E - (1+1i)*tau*Lambda/2);
    u_temp = Lambda*(u(i-1,:))';
    b(2:n_x-1) = u_temp(2:n_x-1) + (f_temp(2:n_x-1))';
    w = A\b;
    u(i,2:n_x-1) = u(i-1,2:n_x-1) + tau*real(w(2:n_x-1)');
end
end
