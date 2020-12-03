function [y, lambda] = PowerIteration(A, eps)
    disp('Power Iteration:');
    n = length(A);
    y = rand(n, 1);
    
    yk = A * y;
    
    y = yk;
    
    lambda = yk(1) / y(1);
    
    inaccuracy = norm(A*y - lambda * y)/norm(y);
    
    iteration = 1;
    
    while inaccuracy > eps
%     while iteration < 100
        yk = A * y;
        lambda = yk(1) / y(1);
        inaccuracy = norm(A*yk - lambda * yk)/norm(yk);
        iteration = iteration + 1;
        y = yk / norm(A * y);
%         disp(inaccuracy);
%         disp(lambda);
    end
    
    disp('Number of iterations');
    disp(iteration);
end
