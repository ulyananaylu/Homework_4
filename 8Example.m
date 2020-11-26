clear; clc; close('all');

%% Inputs:
% specify step
N = 10;
M = 100;
T = 0.1;

[explicitSolution, x, t] = explicitScheme(T, N, M);

exactSolution = zeros(size(explicitSolution));

for i = 1:M+1
    for j = 1:N+1
        exactSolution(i, j) = x(j) + t(i);
    end
end

figure;
surf(x, t, explicitSolution);
title('Explicit scheme');
xlabel('X');
ylabel('T');
zlabel('U');

figure;
surf(x, t, exactSolution);
title('Exact solution');
xlabel('X');
ylabel('T');
zlabel('U');

[implicitSolution, x, t] = implicitScheme(T, N, M);

disp('Explicit deviations');
disp(abs(explicitSolution(M+1,:) - exactSolution(M+1,:))');
disp('Implicit deviations');
disp(abs(implicitSolution(M+1,:) - exactSolution(M+1,:))');

figure;
surf(x, t, implicitSolution);
title('Implicit scheme');
xlabel('X');
ylabel('T');
zlabel('U');
