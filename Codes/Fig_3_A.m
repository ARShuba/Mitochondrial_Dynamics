clear all

%% Parameters for Normal Dynamics
epsilon_normal = 0.01;          % Normal scaled parameter epsilon
epsilon_prime_normal = 0.001;   % Normal scaled parameter epsilon'
q_normal = 0.01;                % Normal parameter q
f_normal = 1;                   % Normal stoichiometric factor

%% Parameters for Cancer-Like Dynamics
epsilon_cancer = 0.005;         % Cancer condition: lower epsilon
epsilon_prime_cancer = 0.0005;  % Cancer condition: lower epsilon'
q_cancer = 0.02;                % Cancer condition: higher q
f_cancer = 0.6;                 % Cancer condition: lower f

% Time span for simulation
tspan = [0 50];  

% Initial conditions for x, y, z
x0 = 0.5;  % Initial fission factor concentration (DRP1)
y0 = 0.2;  % Initial intermediate concentration (MID49/51)
z0 = 1;    % Initial fusion factor concentration (MFN1/2, OPA1)
initial_conditions = [x0, y0, z0];

%% Equations for Normal and Cancer-Like Dynamics
% Normal dynamics
oregonator_eqns_normal = @(t, vars) [
    (1/epsilon_normal) * (q_normal * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
    (1/epsilon_prime_normal) * (-q_normal * vars(2) - vars(1) * vars(2) + f_normal * vars(3));
    vars(1) - vars(3)
];

% Cancer-like dynamics
oregonator_eqns_cancer = @(t, vars) [
    (1/epsilon_cancer) * (q_cancer * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
    (1/epsilon_prime_cancer) * (-q_cancer * vars(2) - vars(1) * vars(2) + f_cancer * vars(3));
    vars(1) - vars(3)
];

[t_normal, vars_normal] = ode45(oregonator_eqns_normal, tspan, initial_conditions);
[t_cancer, vars_cancer] = ode45(oregonator_eqns_cancer, tspan, initial_conditions);

% Extract solutions for normal and cancer-like dynamics
x_normal = vars_normal(:, 1);  % Fission factor (DRP1)
y_normal = vars_normal(:, 2);  % Intermediate (MID49/51)
z_normal = vars_normal(:, 3);  % Fusion factor (MFN1/2, OPA1)

x_cancer = vars_cancer(:, 1);  % Fission factor (DRP1)
y_cancer = vars_cancer(:, 2);  % Intermediate (MID49/51)
z_cancer = vars_cancer(:, 3);  % Fusion factor (MFN1/2, OPA1)

% Log-transformed solutions
x_log_normal = log10(x_normal);
y_log_normal = log10(y_normal);
z_log_normal = log10(z_normal);

x_log_cancer = log10(x_cancer);
y_log_cancer = log10(y_cancer);
z_log_cancer = log10(z_cancer);

%% Figure: Log-Transformed Dynamics
figure;
subplot(3, 1, 1);
plot(t_normal, x_log_normal, 'b', t_cancer, x_log_cancer, 'b--', 'LineWidth', 1.5);
title('Log_{10}(Fission Factor): Normal vs Cancer-Like');
xlabel('Time (\tau)');
ylabel('log_{10}(Concentration)');
legend('Normal', 'Cancer-Like');
grid on;

subplot(3, 1, 2);
plot(t_normal, y_log_normal, 'r', t_cancer, y_log_cancer, 'r--', 'LineWidth', 1.5);
title('Log_{10}(Intermediate): Normal vs Cancer-Like');
xlabel('Time (\tau)');
ylabel('log_{10}(Concentration)');
legend('Normal', 'Cancer-Like');
grid on;

subplot(3, 1, 3);
plot(t_normal, z_log_normal, 'g', t_cancer, z_log_cancer, 'g--', 'LineWidth', 1.5);
title('Log_{10}(Fusion Factor): Normal vs Cancer-Like');
xlabel('Time (\tau)');
ylabel('log_{10}(Concentration)');
legend('Normal', 'Cancer-Like');
grid on;

sgtitle('Log-Transformed Dynamics: Normal vs Cancer-Like');

