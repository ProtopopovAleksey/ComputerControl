
%% Prepare environment
clear all;
close all;
clc;

%% Parse initial condition
fid = fopen('C:\Computer_Controle\input.txt');

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

%fprintf("X: %d\n", X);
%fprintf("N: %d\n", N);

%% Create state matrices
A = rand(size(X,1),size(X,1));
A = [[-0.03333333 0.1 0.01666667 0.2 0.26666667 0.02 ]
[ 0.33333333 0.03333333 -0.33333333 0.33333333 0.06666667 0.33333333]
[-0.33333333 -0.33333333 -0.2 0.33333333 0.1 0.16666667]
[ 0.02333333 -0.03 0.33333333 0.06666667 -0.16666667 -0.26666667]
[ 0.33333333 -0.33333333 0.23333333 0.3 -0.3 -0.1 ]
[ 0.33333333 0.33333333 0.26666667 0.16666667 0.13333333 0.06666667]];

B = randn(size(X,1),size(U,1));
E = randn(size(X,1),size(N,1));
C = randn(size(Y,1),size(X,1));
D = randn(size(Y,1),size(U,1));
F = randn(size(Y,1),size(N,1));

%% Determine stability
self_values = eig(A);
stability = 1;
for i=1:size(self_values,1) 
    if (abs(self_values(i))>1)
       stability = 0; 
    elseif (abs(self_values(i))==1)
       stability = 2; 
    end
end
%fprintf("self_values: %d\n", self_values);
if (stability == 0)
    fprintf("System instable\n"); 
elseif (stability == 1)
	fprintf("System stable\n");
elseif (stability == 2)
	fprintf("System at the stability boundary\n");
end

%% Modeling system
X_log = zeros(size(X,1),T);
Y_log = zeros(size(Y,1),T);
% Modelling
for i = 1:T
    X = A*X + B*U + E*N;
    Y = C*X + D*U + F*N;
    % Save history
    X_log(:,i) = X;
    Y_log(:,i) = Y;
end

% Plot graph
% X states
figure(1);
plot(X_log');
grid on;
title('X');

% Y states
figure(2);
plot(Y_log');
grid on;
title('Y');

%% Determine step parameters

% Determine overshoot
overshoot = (max(Y_log) - Y_log(T))/Y_log(T)*100;
%fprintf("final_value: %d\n", Y_log(T));
%fprintf("max(Y_log): %d\n", max(Y_log));
fprintf("overshoot: %d\n", overshoot);

% Determine step response time
step_time = 0;
for i=1:length(Y_log)
    if abs(Y_log(i)-Y_log(T))/abs(Y_log(T)) > 0.05
        step_time=i;
    end    
end
fprintf("step_time: %d\n", step_time);

%% Filtering

%% Modeling filtering system
X_log_f = zeros(size(X,1),T);
Y_log_f = zeros(size(Y,1),T);

alpha = 0.05;
% Modelling
for i = 1:T
        % First step
    if(i == 1)
        X = A*X + B*U + E*N;
        X_f = X;
    else
        % Filter
        X_f = X_f + alpha*(X - X_f);
        X = A*X_f + B*U + E*N;
    end
    Y = C*X + D*U + F*N;
    % Save history
    X_log_f(:,i) = X;
    Y_log_f(:,i) = Y;
end

% Plot graph
% X states
figure(3);
plot(X_log_f');
grid on;
title('X_f');

% Y states
figure(4);
plot(Y_log_f');
grid on;
title('Y_f');


