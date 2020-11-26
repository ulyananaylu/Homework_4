function [y, lambda] = ScalarProduct(A, eps)
    disp('Scalar Product:');
    n = length(A);
    y = rand(n, 1);
    
    lambda = (y' * A * y)/(y' * y);
    
    inaccuracy = norm(A*y - lambda * y)/norm(y);
    
    iteration = 1;
    
    while inaccuracy > eps 
        y = (A * y) / norm(A * y);
        lambda = (y' * A * y)/(y' * y);
        inaccuracy = norm(A*y - lambda * y)/norm(y);
        iteration = iteration + 1;
    end
    
    disp('Number of iterations');
    disp(iteration);
end
