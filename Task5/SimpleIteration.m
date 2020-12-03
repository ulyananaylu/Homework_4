function [x] = SimpleIteration(A, b, eps)
    disp('Simple iterations:');
    D = diag(diag(A));
    
    H = eye(size(A)) - inv(D) * A;
    
    x0 = zeros(length(b), 1);
    
    if (max(abs(eig(H))) < 1)
        g = inv(D) * b;
    
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
