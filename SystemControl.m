%% Очистка
clear all;
close all;
clc;

%% Чтение из файла
fid = fopen('C:\Computer_Control\input.txt');

tline = fgetl(fid);
X = str2num(tline);
tline = fgetl(fid);
U = str2num(tline);
tline = fgetl(fid);
N = str2num(tline);
tline = fgetl(fid);
Y = str2num(tline);
tline = fgetl(fid);
T = str2num(tline);

X= X';
U= U';
N= N';
Y= Y';

fclose(fid);
%% Создание матриц
%A = rand(size(X,1),size(X,1));
A = [[-0.03 0.1 0.016 0.2 0.26 0.02 ]
[ 0.3 0.03 -0.3 0.3 0.06 0.3]
[-0.3 -0.3 -0.2 0.3 0.1 0.16]
[ 0.02 -0.03 0.3 0.06 -0.16 -0.26]
[ 0.3 -0.33 0.23 0.3 -0.3 -0.1 ]
[ 0.3 0.3 0.26 0.16 0.13 0.06]];

B = randn(size(X,1),size(U,1));
E = randn(size(X,1),size(N,1));
C = randn(size(Y,1),size(X,1));
D = randn(size(Y,1),size(U,1));
F = randn(size(Y,1),size(N,1));


%% Определение устойчивости
self_values = eig(A);
stability = 1;
for i=1:size(self_values,1) 
    if (abs(self_values(i))>1)
       stability = 0; 
    elseif (abs(self_values(i))==1)
       stability = 2; 
    end
end

if (stability == 0)
    fprintf("Система неустойчива\n"); 
elseif (stability == 1)
	fprintf("Система устойчива\n");
elseif (stability == 2)
	fprintf("Система на границе устойчивости\n");
end

%% Моделирование системы + экспоненциальный фильтр
X_log = zeros(size(X,1),T);
Y_log = zeros(size(Y,1),T);

alpha = 0.2;

for i = 1:T
    if(i == 1)
        X = A*X + B*U + E*N;
        Y = C*X + D*U + F*N;
        Y_f = Y;
    else
        X = A*X + B*U + E*N;
        Y = C*X + D*U + F*N;
        Y_f(i) = alpha*Y + (1-alpha)*Y_f(i-1);
    end

    X_log(:,i) = X;
    Y_log(:,i) = Y;
end

dlmwrite('C:\Computer_Control\output.txt',  Y_log);

figure(1);
plot(X_log');
grid on;
title('X');


figure(2);
plot(Y_log');
grid on;
title('Y');

figure(3);
plot(Y_f);
grid on;
title('Y_f');

%% Определение параметров переходного процесса
fprintf("Параметры переходного процесса:\n"); 
% Перерегулирование
overshoot = (max(Y_log) - Y_log(T))/Y_log(T)*100;
fprintf("overshoot: %d\n", overshoot);

% Время переходного процесса
step_time = 0;
for i=1:length(Y_log)
    if abs(Y_log(i)-Y_log(T))/abs(Y_log(T)) > 0.05
        step_time=i;
    end    
end
fprintf("step_time: %d\n", step_time);

%% Определение параметров переходного процесса отфильтрованной системы

fprintf("Параметры переходного процесса с фильтром:\n"); 
% Перерегулирование
overshoot = (max(Y_f) - Y_f(T))/Y_f(T)*100;
fprintf("overshoot_f: %d\n", overshoot);

% Время переходного процесса
step_time = 0;
for i=1:length(Y_log)
    if abs(Y_f(i)-Y_f(T))/abs(Y_f(T)) > 0.05
        step_time=i;
    end    
end
fprintf("step_time_f: %d\n", step_time);


