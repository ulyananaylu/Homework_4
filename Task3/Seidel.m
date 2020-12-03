function [x] = Seidel(A, b, eps)
    disp('Seidel:');
    D = diag(diag(A));
    
    L = tril(A, -1);
    R = triu(A, 1);
    
    H = zeros(size(A)) - inv(D + L) * R;
    g = inv(D + L) * b;
    
    
    x0 = zeros(length(b), 1);
    
    if (max(abs(eig(H))) < 1)
        x = H*x0 + g;
    
        iteration = 1;
    
        while norm(x - x0) > eps
            iteration = iteration + 1;
        
            c = x;
        
            x = H * x + g;
        
            x0 = c;
        end
    
        disp('Number of iterations:');
        disp(iteration);
    else
        disp('Error: H matrix spectral radius is grater than 1');
        x = x0;
    end
end
