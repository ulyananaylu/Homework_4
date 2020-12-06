function [y, t, d] = Richardson(method, t0, y0, tn, A)
    if method == "Adams3"
        h = 1/(max(abs(eig(A))));
        order = 2;
    elseif method == "RungeKutta"
        h = 2.78/(max(abs(eig(A))));
        order = 4;
    elseif method == "BackwardEuler"
        h = 1e-3;
        order = 2;
    end
    epsilon = 1e-5;
    [t1, v1] = feval(method,t0, y0, h, tn, A);
    d = ones(1,rank(A));
    while max(abs(d))>epsilon
        h = h/2;
        [t2, v2] = feval(method,t0, y0, h, tn, A);
        dd = zeros(3,rank(A));
        dd(1,:) = (v2(:,1) - v1(:,1))/(2^order - 1);
        m = 0;
        for j = 2:length(t1)
            dd(3,:) = (v2(:,2*j-1) - v1(:,j))/(2*order-1);
            dd(2,:) = (dd(1,:) + dd(3,:))/2;
            for i = 1:3
                if max(abs(dd(i,:)))>m
                    m = max(abs(dd(i,:)));
                    d = dd(i,:);
                end
            end
            dd(1,:) = dd(3,:);
        end
        fprintf(1, '%1u %2.5g \n', length(t2),m);
        v1 = v2; t1 = t2;
    end
    y = v2;
    t = t2;
end
