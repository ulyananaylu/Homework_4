function [x,t, u] = tochn(a,T,tau,h,f)
x = 0:h:a;
t = 0:tau:T;
n_x = length(x);
n_t = length(t);
u_resh = zeros(n_t,n_x);
for i = 1:n_t
    for j = 1:n_x
        u_resh(i,j) = f(x(j),t(i));
    end
end
u = u_resh;
end