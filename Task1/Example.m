clear; clc; close('all');
A = [
        -125    123.55;
        123.55  -123
    ];
y0 = [
        1; 
        1
     ];
 
t0 = 0;
h = 0.001; % try: h = 0.01, 0.001
tn = 0.5; % equal to: t0 + h*n, with n the number of steps

[t, y] = Euler(t0, y0, h, tn, A);

figure;
plot(t, y(1,:), 'b');
title('Euler method');
hold('on');
plot(t, y(2,:), 'r');
legend('Y1', 'Y2');

[t, y] = BackwardEuler(t0, y0, h, tn, A, true);

figure;
plot(t, y(1,:), 'b');
title('Backward Euler');
hold('on');
plot(t, y(2,:), 'r');
legend('Y1', 'Y2');

[t, y] = Adams3(t0, y0, h, tn, A);

figure;
plot(t, y(1,:), 'b');
title('Third order Adams method');
hold('on');
plot(t, y(2,:), 'r');
legend('Y1', 'Y2');

[t, y] = RungeKutta(t0, y0, h, tn, A);

figure;
plot(t, y(1,:), 'b');
title('Runge-Kutta');
hold('on');
plot(t, y(2,:), 'r');
legend('Y1', 'Y2');
