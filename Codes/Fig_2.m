clear all
%%
% MATLAB Code for Mitochondrial Dynamics Figures

% Parameters
epsilon = 0.01;          % Scaled parameter epsilon
epsilon_prime = 0.001;    % Scaled parameter epsilon'
q = 0.01;               % Parameter q
f = 1;                    % Stoichiometric factor

% % Parameters
% epsilon = 0.01;          % Scaled parameter epsilon
% epsilon_prime = 0.001;   % Scaled parameter epsilon'
% q = 0.01;                % Parameter q
% f = 1;                   % Stoichiometric factor

% Time span for simulation
tspan = [0 50]; 

% Initial conditions for x, y, z
x0 = 0.5;  % Initial fission factor concentration (DRP1)
y0 = 0.2;  % Initial intermediate concentration (MID49/51)
z0 = 1;    % Initial fusion factor concentration (MFN1/2, OPA1)
initial_conditions = [x0, y0, z0];

% scaled Field-Noyes equations
oregonator_eqns = @(t, vars) [
    (1/epsilon) * (q * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));  % dx/dτ
    (1/epsilon_prime) * (-q * vars(2) - vars(1) * vars(2) + f * vars(3));      % dy/dτ
    vars(1) - vars(3)                                                         % dz/dτ
];

[t, vars] = ode45(oregonator_eqns, tspan, initial_conditions);

% Extract solutions
x = vars(:, 1);  % Fission factor (DRP1)
y = vars(:, 2);  % Intermediate (MID49/51)
z = vars(:, 3);  % Fusion factor (MFN1/2, OPA1)

% log transformed solutions
x_log = log10(x);
y_log = log10(y);
z_log = log10(z);

%% Figure 2A-B: Phase Plane Plot (Subplots)
figure;
subplot(2, 1, 1);
plot(x, y, 'k', 'LineWidth', 1.5);
title('Phase Plane: Fission Factor (x) vs Intermediate (y)');
xlabel('Fission Factor (x)');
ylabel('Intermediate (y)');
grid on;

subplot(2, 1, 2);
plot(x, z, 'm', 'LineWidth', 1.5);
title('Phase Plane: Fission Factor (x) vs Fusion Factor (z)');
xlabel('Fission Factor (x)');
ylabel('Fusion Factor (z)');
grid on;

sgtitle('Mitochondrial Dynamics: Phase Planes');

%% Figure 2C: Parameter Sweep (Subplots)
% Parameter sweeps for epsilon, epsilon_prime, and f
epsilon_values = [0.005, 0.01, 0.02, 0.03];
epsilon_prime_values = [0.0005, 0.001, 0.002];
f_values = [0.8, 1, 1.2];

figure;

% Sweep epsilon
subplot(3, 1, 1);
hold on;
for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    oregonator_eqns = @(t, vars) [
        (1/epsilon) * (q * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
        (1/epsilon_prime) * (-q * vars(2) - vars(1) * vars(2) + f * vars(3));
        vars(1) - vars(3)
    ];
    [t, vars] = ode45(oregonator_eqns, tspan, initial_conditions);
    plot(t, vars(:, 1), 'LineWidth', 1.5);
end
title('Sweep: \epsilon');
xlabel('Time (\tau)');
ylabel('Fission Factor (x)');
legend('\epsilon = 0.005', '\epsilon = 0.01', '\epsilon = 0.02', '\epsilon = 0.03');
grid on;

% Sweep epsilon_prime
subplot(3, 1, 2);
hold on;
for i = 1:length(epsilon_prime_values)
    epsilon_prime = epsilon_prime_values(i);
    oregonator_eqns = @(t, vars) [
        (1/epsilon) * (q * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
        (1/epsilon_prime) * (-q * vars(2) - vars(1) * vars(2) + f * vars(3));
        vars(1) - vars(3)
    ];
    [t, vars] = ode45(oregonator_eqns, tspan, initial_conditions);
    plot(t, vars(:, 1), 'LineWidth', 1.5);
end
title('Sweep: \epsilon''');
xlabel('Time (\tau)');
ylabel('Fission Factor (x)');
legend('\epsilon'' = 0.0005', '\epsilon'' = 0.001', '\epsilon'' = 0.002');
grid on;

% Sweep f
subplot(3, 1, 3);
hold on;
for i = 1:length(f_values)
    f = f_values(i);
    oregonator_eqns = @(t, vars) [
        (1/epsilon) * (q * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));
        (1/epsilon_prime) * (-q * vars(2) - vars(1) * vars(2) + f * vars(3));
        vars(1) - vars(3)
    ];
    [t, vars] = ode45(oregonator_eqns, tspan, initial_conditions);
    plot(t, vars(:, 1), 'LineWidth', 1.5);
end
title('Sweep: Stoichiometric Factor (f)');
xlabel('Time (\tau)');
ylabel('Fission Factor (x)');
legend('f = 0.8', 'f = 1', 'f = 1.2');
grid on;

sgtitle('Parameter Sweeps for \epsilon, \epsilon'', and f');
