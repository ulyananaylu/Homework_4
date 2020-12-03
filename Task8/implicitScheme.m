function [u, x, t] = implicitScheme(T, N, M)
    h = 1 / N;
    tau = T / M;
    u = zeros(M + 1, N + 1);
    
    x = (0: h: 1);
    t = (0: tau: T)';
    
    for i = 1:length(x)
        u(1,i) = phi(x(i));
    end
    
    A = zeros(1, N + 1);
    B = zeros(1, N + 1);
    C = zeros(1, N + 1);
    G = zeros(1, N + 1);
    
    for k = 2 : M + 1 
        [alpha1, alpha2, alpha] = leftBoundary(t(k));
        [beta1, beta2, beta] = rightBoundary(t(k));
        
        B(1) =  -alpha1 - alpha2/h;
        C(1) = -alpha2/h;
        G(1) = alpha;
        
        B(N + 1) = - beta1 - beta2/h;
        A(N + 1) = -beta2/h;
        G(N + 1) = beta;
        
        for i = 2 : N 
            C(i) = a(x(1), t(1)) / h^2 + b(x(1), t(1)) / (2 * h);
            B(i) =  2 * a(x(1), t(1)) / h^2 + 1/tau + c(x(1), t(1));
            A(i) = a(x(1), t(1)) / h^2 - b(x(1), t(1)) / (2 * h);
            G(i) = - u(k-1, i)/tau - f(x(i-1), t(k-1));
        end
        
        u(k,:) = tridiagonalMatrix(A, B, C, G);
    end
end
