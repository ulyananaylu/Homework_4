function [u, x, t] = explicitScheme(T, N, M)
    h = 1 / N;
    tau = T / M;
    u = zeros(M + 1, N + 1);
    
    x = (0: h: 1);
    t = (0: tau: T)';
    
    A = zeros(length(t), length(x));
    for i = 1:length(t)
        for j = 1:length(x)
            A(i,j) = a(x(j),t(i));
        end
    end
    
    if (max(max(A)) * tau / h^2 <= 1/2) 
        for i = 1:length(x)
            u(1,i) = phi(x(i));
        end

        for k = 2:M + 1
            for i = 2:N
                u(k, i) = u(k-1, i) + tau * (Lh(x(i), t(k - 1), u(k -1, i+1), u(k - 1,i), u(k - 1, i-1), h) + f(x(i), t(k - 1)));
            end
            [alpha1, alpha2, alpha] = leftBoundary(t(k));
            [beta1, beta2, beta] = rightBoundary(t(k));
            
            u(k, 1) = (alpha + alpha2 * (4 * u(k, 2) - u(k, 3)) / (2*h)) / (alpha1 + (3 * alpha2) / (2 * h));
            u(k, N + 1) = (beta + beta2 * (4 * u(k, N) - u(k, N -1)) / (2 * h)) / (beta1 + (3 * beta2) / (2 * h));
        end
    end
end
