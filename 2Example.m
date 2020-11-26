clear; clc; close('all');

%% Inputs:
%% Matrix
A = [
        7.35272     0.88255     -2.170052;
        0.88255     5.58351     0.528167;
        -2.170052    0.5528167    4.430329
     ];
 
 %% Right vector
 y = [
        1;
        0;
        0
     ];
%% Solve using LU defactorization


format long;
[L, U] = LU(A);

disp(L);
disp(U);
disp(L*U);

x = LUsystem(A, y);
disp('Answer:');
disp(x);

%% Solve using square root method
x = Gsystem(A, y);

disp('Answer:');
disp(x);

x = A \ y;
disp(x);
