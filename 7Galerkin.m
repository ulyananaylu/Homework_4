function [y] = Galerkin(n, a, b)    
    [phi, dphi, ddphi] = myphi_dphi_d2kk_good_anal(1, n);
    
    A = zeros(n);
    
    f = zeros(n, 1);
    
    for i = 1:n
        for j = 1:n
            A(i,j) = int((P * ddphi(j) + Q * dphi(j) + R * phi(j)) * phi(i),a,b) ;
        end
        f(i) = int(F * phi(i), a, b);
    end    
    
    c = A \ f;
    
    y = phi' * c;
end
