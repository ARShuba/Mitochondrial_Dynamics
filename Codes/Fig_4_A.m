clear all
%% Parameters

% Baseline Parameters (Disease State)
epsilon_disease = 0.005;         % Disease condition: lower epsilon
epsilon_prime_disease = 0.0005;  % Disease condition: lower epsilon'
q_disease = 0.02;                % Disease condition: higher q
f_disease = 0.6;                 % Disease condition: lower stoichiometric factor

% Therapeutic Intervention Parameters
f_therapeutic = 1.2;  % Enhanced fusion promotion (increased z dynamics)

% Time span for simulation
tspan = [0 50];  

% Initial conditions for x, y, z
x0 = 0.5;  % Initial fission factor concentration (DRP1)
y0 = 0.2;  % Initial intermediate concentration (MID49/51)
z0 = 1;    % Initial fusion factor concentration (MFN1/2, OPA1)
initial_conditions = [x0, y0, z0];

%% Equations for Disease and Therapeutic States

% Disease state dynamics
oregonator_eqns_disease = @(t, vars) [
    (1/epsilon_disease) * (q_disease * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
    (1/epsilon_prime_disease) * (-q_disease * vars(2) - vars(1) * vars(2) + f_disease * vars(3));
    vars(1) - vars(3)
];

% Therapeutic intervention dynamics
oregonator_eqns_therapeutic = @(t, vars) [
    (1/epsilon_disease) * (q_disease * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
    (1/epsilon_prime_disease) * (-q_disease * vars(2) - vars(1) * vars(2) + f_therapeutic * vars(3));
    vars(1) - vars(3)
];


[t_disease, vars_disease] = ode45(oregonator_eqns_disease, tspan, initial_conditions);
[t_therapeutic, vars_therapeutic] = ode45(oregonator_eqns_therapeutic, tspan, initial_conditions);

% Extract solutions for disease and therapeutic states
x_disease = vars_disease(:, 1);  % Fission factor (DRP1)
y_disease = vars_disease(:, 2);  % Intermediate (MID49/51)
z_disease = vars_disease(:, 3);  % Fusion factor (MFN1/2, OPA1)

x_therapeutic = vars_therapeutic(:, 1);  % Fission factor (DRP1)
y_therapeutic = vars_therapeutic(:, 2);  % Intermediate (MID49/51)
z_therapeutic = vars_therapeutic(:, 3);  % Fusion factor (MFN1/2, OPA1)

% log transformed solutions
x_disease_log = log10(x_disease);
y_disease_log = log10(y_disease);
z_disease_log = log10(z_disease);

x_therapeutic_log = log10(x_therapeutic);
y_therapeutic_log = log10(y_therapeutic);
z_therapeutic_log = log10(z_therapeutic);

%% Figure: Log-transformed dynamics (Comparative Time Series)
figure;
subplot(3, 1, 1);
plot(t_disease, x_disease_log, 'b', t_therapeutic, x_therapeutic_log, 'b--', 'LineWidth', 1.5);
title('Log_{10}(Fission Factor, DRP1): Disease vs Therapeutic');
xlabel('Time (\tau)');
ylabel('Log_{10}(x)');
legend('Disease State', 'Therapeutic Intervention');
grid on;

subplot(3, 1, 2);
plot(t_disease, y_disease_log, 'r', t_therapeutic, y_therapeutic_log, 'r--', 'LineWidth', 1.5);
title('Log_{10}(Intermediate, MID49/51): Disease vs Therapeutic');
xlabel('Time (\tau)');
ylabel('Log_{10}(y)');
legend('Disease State', 'Therapeutic Intervention');
grid on;

subplot(3, 1, 3);
plot(t_disease, z_disease_log, 'g', t_therapeutic, z_therapeutic_log, 'g--', 'LineWidth', 1.5);
title('Log_{10}(Fusion Factor, MFN1/2, OPA1): Disease vs Therapeutic');
xlabel('Time (\tau)');
ylabel('Log_{10}(z)');
legend('Disease State', 'Therapeutic Intervention');
grid on;

sgtitle('Log-transformed Dynamics: Disease vs Therapeutic Intervention');
