% Лабораторная работа №1
% Решение задачи о назначениях венгерским методом

clear; clc;

% Параметры
debug = true;
maximization = false;

% Исходная матрица стоимостей
C = [3, 5, 2, 4, 8;
     10, 10, 4, 3, 6;
     5, 6, 9, 8, 3;
     6, 2, 5, 8, 4;
     5, 4, 8, 9, 3];

[X_opt, f_opt] = hungarian(C, debug, maximization);

disp("Оптимальная матрица назначений X_opt = ");
disp(X_opt);

fprintf('f_opt: %g\n', f_opt);

