clear all

%%
% Parameters
epsilon = 3.6e-2;          % Scaled parameter epsilon
epsilon_prime = 1.2e-4;    % Scaled parameter epsilon'
q = 2.4e-4;               % Parameter q
f = 1;                    % Stoichiometric factor

% Time span for simulation
tspan = [0 50];  

% Initial conditions for x, y, z
x0 = 1;  % Initial activator concentration
y0 = 1;  % Initial intermediate concentration
z0 = 1;  % Initial inhibitor concentration
initial_conditions = [x0, y0, z0];

% scaled Field-Noyes equations
oregonator_eqns = @(t, vars) [
    (1/epsilon) * (q * vars(2) - vars(1) * vars(2) + vars(1) * (1 - vars(1)));  % dx/dτ
    (1/epsilon_prime) * (-q * vars(2) - vars(1) * vars(2) + f * vars(3));      % dy/dτ
    vars(1) - vars(3)                                                         % dz/dτ
];

[t, vars] = ode45(oregonator_eqns, tspan, initial_conditions);

% Extract solutions
x = vars(:, 1);  % Activator
y = vars(:, 2);  % Intermediate
z = vars(:, 3);  % Inhibitor

% Log transformed solution
x_log = log10(x);
y_log = log10(y);
z_log = log10(z);

%% Figure 1B: Log-Transformed Oscillations (Subplots)
figure;
subplot(3, 1, 1);
plot(t, x_log, 'b', 'LineWidth', 1.5);
title('Log_{10}(Activator) Concentration');
xlabel('Time (\tau)');
ylabel('log_{10}(x)');
grid on;

subplot(3, 1, 2);
plot(t, y_log, 'r', 'LineWidth', 1.5);
title('Log_{10}(Intermediate) Concentration');
xlabel('Time (\tau)');
ylabel('log_{10}(y)');
grid on;

subplot(3, 1, 3);
plot(t, z_log, 'g', 'LineWidth', 1.5);
title('Log_{10}(Inhibitor) Concentration');
xlabel('Time (\tau)');
ylabel('log_{10}(z)');
grid on;

sgtitle('Log-Transformed Oscillations');