clear all
%%

% Parameters
epsilon = 0.01;          % Scaled parameter epsilon
epsilon_prime = 0.001;    % Scaled parameter epsilon'
q = 0.01;               % Parameter q
f = 1;                    % Stoichiometric factor


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


%% Figure 1D: Log-Transformed Oscillations (Subplots)
figure;
subplot(3, 1, 1);
plot(t, x_log, 'b', 'LineWidth', 1.5);
title('Log_{10}(Fission Factor) (DRP1)');
xlabel('Time (\tau)');
ylabel('log_{10}(x)');
grid on;

subplot(3, 1, 2);
plot(t, y_log, 'r', 'LineWidth', 1.5);
title('Log_{10}(Intermediate) (MID49/51)');
xlabel('Time (\tau)');
ylabel('log_{10}(y)');
grid on;

subplot(3, 1, 3);
plot(t, z_log, 'g', 'LineWidth', 1.5);
title('Log_{10}(Fusion Factor) (MFN1/2, OPA1)');
xlabel('Time (\tau)');
ylabel('log_{10}(z)');
grid on;

sgtitle('Mitochondrial Dynamics: Log-Transformed Oscillations');