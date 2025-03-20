% Furuta Pendulum Simulation (Nonlinear Dynamics with Inputs)

clear; clc; close all;
prm.eps = 0.08; % deviation from initial condition
prm.dist = +0.1; %intensity of the disturbance torque
%% Parameters (Original System)
 % System physical parameters
prm.L1 = 0.278;    % Length parameter 1 [m]
prm.L2 = 0.300;    % Length parameter 2 [m]
prm.l1 = 0.150;    % Center of mass distance for first link [m]
prm.l2 = 0.148;    % Center of mass distance for second link [m]
prm.m1 = 0.300;    % Mass of first link [kg]
prm.m2 = 0.075;    % Mass of second link [kg]
prm.J1 = 2.48e-2;  % Inertia of first link [kg*m^2]
prm.J2 = 3.86e-3;  % Inertia of second link [kg*m^2]
prm.b1 = 1.00e-4;  % Damping coefficient for first link
prm.b2 = 2.80e-4;  % Damping coefficient for second link
prm.g  = 9.81;     % Gravitational acceleration [m/s^2]
prm.Km = 0.090;    % Motor constant
prm.Lm = 0.005;    % Motor inductance [H]
prm.Rm = 7.80;     % Motor resistance [Ohm]

prm.J0_hat = prm.J1 + prm.m1*prm.l1^2 + prm.m2*prm.L1^2;
prm.J2_hat = prm.J2 + prm.m2*prm.l2^2;

%% Define Control Inputs (Torque Inputs)
V = @(t) 0;  % Voltage applied to the motor (V)
u2 = @(t) 0;  % External disturbance torque on Arm 2 (Nm)

%% Equations of Motion with Inputs
pendulum_dynamics = @(t, y) [
    y(3);
    y(4);
    (- prm.J2_hat*prm.b1*y(3) ...
     - prm.m2*prm.L1*prm.l2*cos(y(2))*prm.b2*y(4) ...
     - prm.J2_hat^2*sin(2*y(2))*y(3)*y(4) ...
     + 0.5*prm.J2_hat*prm.m2*prm.L1*prm.l2*cos(y(2))*sin(2*y(2))*y(3)^2 ...
     - prm.J2_hat*prm.m2*prm.L1*prm.l2*sin(y(2))*y(4)^2 ...
     + prm.J2_hat*prm.Km*y(5) ...
     + prm.m2*prm.L1*prm.l2*cos(y(2))*u2(t) ...
     + 0.5*prm.m2^2*prm.l2^2*prm.L1*sin(2*y(2))*prm.g) ...
     /(prm.J0_hat*prm.J2_hat + prm.J2_hat^2*sin(y(2))^2 - prm.m2^2*prm.L1^2*prm.l2^2*cos(y(2))^2);

    (- prm.m2*prm.L1*prm.l2*cos(y(2))*prm.b1*y(3) ...
     - prm.b2*(prm.J0_hat+prm.J2_hat*sin(y(2))^2)*y(4) ...
     - prm.m2*prm.L1*prm.l2*prm.J2_hat*cos(y(2))*sin(2*y(2))*y(3)*y(4) ...
     - 0.5*sin(2*y(2))*(prm.J0_hat*prm.J2_hat + prm.J2_hat^2*sin(y(2))^2)*y(3)^2 ...
     - 0.5*prm.m2^2*prm.L1^2*prm.l2^2*sin(2*y(2))*y(4)^2 ...
     + prm.m2*prm.L1*prm.l2*cos(y(2))*prm.Km*y(5) ...
     + (prm.J0_hat+prm.J2_hat*sin(y(2))^2)*u2(t) ...
     + prm.m2*prm.l2*sin(y(2))*(prm.J0_hat+prm.J2_hat*sin(y(2))^2)*prm.g)...
     /(prm.J0_hat*prm.J2_hat + prm.J2_hat^2*sin(y(2))^2 - prm.m2^2*prm.L1^2*prm.l2^2*cos(y(2))^2);

    (V(t) - prm.Rm*y(5) - prm.Km*y(3)) / prm.Lm;
];

%% Simulation Settings
Tspan = [0 100];
Y0 = [0; pi; 0; 0; 0]; % Initial conditions: [theta1, theta2, theta1_dot, theta2_dot]

%% Solve ODE using a more stable solver
[t, Y] = ode15s(pendulum_dynamics, Tspan, Y0);

%% Plot Results
figure(1);
plot(t, Y(:,1), 'r', t, +Y(:,2), 'b');
legend('\theta_1 (Arm Rotation)', '\theta_2 (Pendulum Angle)');
xlabel('Time (s)'); ylabel('Angle (rad)'); title('Furuta Pendulum Simulation with Input Torque');

grid on;
%% LINEARIZATION OF THE SYSTEM 
for i =1:1
%% Define symbolic variables for linearization
syms y1 y2 y3 y4 y5 V u2 real
y = [y1; y2; y3; y4; y5];  % State vector [theta1, theta2, theta1_dot, theta2_dot]
u = [V; u2];  % Input vector [u1, u2]

%% Equations of Motion (Nonlinear)
dy1 = y3;
dy2 = y4;
dy3 = (- prm.J2_hat*prm.b1*y3 ...
     - prm.m2*prm.L1*prm.l2*cos(y2)*prm.b2*y4 ...
     - prm.J2_hat^2*sin(2*y2)*y3*y4 ...
     + 0.5*prm.J2_hat*prm.m2*prm.L1*prm.l2*cos(y2)*sin(2*y2)*y3^2 ...
     - prm.J2_hat*prm.m2*prm.L1*prm.l2*sin(y2)*y4^2 ...
     + prm.J2_hat*prm.Km*y5 ...
     + prm.m2*prm.L1*prm.l2*cos(y2)*u2 ...
     + 0.5*prm.m2^2*prm.l2^2*prm.L1*sin(2*y2)*prm.g) ...
     /(prm.J0_hat*prm.J2_hat + prm.J2_hat^2*sin(y2)^2 - prm.m2^2*prm.L1^2*prm.l2^2*cos(y2)^2);

dy4 = (- prm.m2*prm.L1*prm.l2*cos(y2)*prm.b1*y3 ...
     - prm.b2*(prm.J0_hat+prm.J2_hat*sin(y2)^2)*y4 ...
     - prm.m2*prm.L1*prm.l2*prm.J2_hat*cos(y2)*sin(2*y2)*y3*y4 ...
     - 0.5*sin(2*y2)*(prm.J0_hat*prm.J2_hat + prm.J2_hat^2*sin(y2)^2)*y3^2 ...
     - 0.5*prm.m2^2*prm.L1^2*prm.l2^2*sin(2*y2)*y4^2 ...
     + prm.m2*prm.L1*prm.l2*cos(y2)*prm.Km*y5 ...
     + (prm.J0_hat+prm.J2_hat*sin(y2)^2)*u2 ...
     + prm.m2*prm.l2*sin(y2)*(prm.J0_hat+prm.J2_hat*sin(y2)^2)*prm.g)...
     /(prm.J0_hat*prm.J2_hat + prm.J2_hat^2*sin(y2)^2 - prm.m2^2*prm.L1^2*prm.l2^2*cos(y2)^2);

dy5 = (V - prm.Rm*y5 - prm.Km*y3) / prm.Lm;

% System of equations
f = [dy1; dy2; dy3; dy4; dy5];

%% Linearize the system around equilibrium point
% The equilibrium point is y1=0, y2=pi, y3=0, y4=0, y5 = 0;
y_eq = [0; 0; 0; 0; 0];  % Equilibrium point [theta1, theta2, theta1_dot, theta2_dot, i_motor]
u_eq = [0; 0];  % Equilibrium input [u1, u2]

% Jacobian of the system with respect to state variables and inputs
A = jacobian(f, y);  % Jacobian with respect to the state vector y
B = jacobian(f, u);  % Jacobian with respect to the input vector u

% Evaluate the Jacobians at the equilibrium point
A_eq = double(subs(A, [y; u], [y_eq; u_eq]));  % Evaluate A at equilibrium
B_eq = double(subs(B, [y; u], [y_eq; u_eq]));  % Evaluate B at equilibrium

% Display the linearized system matrices
disp('Linearized A Matrix:');
disp(A_eq);
disp('Linearized B Matrix:');
disp(B_eq);
end
%% DISCRETIZATION

Ts = 0.001; % Sampling time

% Define the continuous-time system
sys_cont = ss(A_eq, B_eq, eye(5), zeros(5, 2));

% Discretization using zero-order hold (ZOH)
sys_disc = c2d(sys_cont, Ts, 'zoh');

% Extract discrete-time system matrices
A_d = sys_disc.A;
B_d = sys_disc.B;

% Separate control and disturbance input matrices
B_c = B_d(:,1);  % First input is the controllable one
B_distr = B_d(:,2);  % Second input acts as a disturbance

% Reference state (we want x_2 to reach pi)

% Define LQR weighting matrices
Q = [ 1 0 0 0 0;
      0 1 0 0 0;
      0 0 0.01 0 0;
      0 0 0 0.01 0;
      0 0 0 0 10];  % Weighting on states (identity matrix)
R = 100;             % Weighting on control input (scalar, as there's only one control input)

% Solve the discrete-time Riccati equation
[K, P, S] = dlqr(A_d, B_c(:,1), Q, R); % Optimal feedback gain

% Simulation parameters
T_sim = 30;
N = T_sim/Ts; % Number of time steps
x_lin = zeros(5, N); % Initialize state trajectory
x_lin(:,1) = [0; 6*pi/180; 0; 0; 0]; % Initial conditions

% Disturbance affecting the second input
d = zeros(1, N);
d(1, N/2:N/2+1) = 0;  % Introduce an impulse disturbance at half simulation

% Closed-loop system simulation with disturbance
for k = 1:N-1
    V = max(min(K*x_lin(:,k),10),-10);
    x_lin(:,k+1) = A_d *x_lin(:,k) - B_c* V + B_distr * d(:,k); % State evolution
end

% Plot the states
figure(2);
subplot(2,2,1);
stairs(0:Ts:T_sim-Ts, x_lin(1,:));
title('State x_1 (\theta_1)');
grid on;

subplot(2,2,2);
stairs(0:Ts:T_sim-Ts, x_lin(2,:));
title('State x_2 (\theta_2)');
grid on;
yline(pi, '--r', 'Reference \pi');  % Reference line for theta2

subplot(2,2,3);
stairs(0:Ts:T_sim-Ts, x_lin(3,:));
title('State x_3 (\theta_1 dot)');
grid on;

subplot(2,2,4);
stairs(0:Ts:T_sim-Ts, x_lin(4,:));
title('State x_4 ({\theta}_2 dot)');
grid on;

xlabel('Time (s)');
sgtitle('State Evolution with LQR Control Linearized system');

% Graphic simulation
% 
% figure;
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% view(3);
% xlim([-0.5, 0.5]); ylim([-0.5, 0.5]); zlim([-0.2, 0.5]);
% title('Simulazione Grafica del Furuta Pendulum con Controllore LQR');
% 
% % Parametri per la visualizzazione
% arm_length = L1;         % Lunghezza dell'arto orizzontale
% pendulum_length = L2;    % Lunghezza del pendolo
% 
% % Inizializza oggetti grafici (con posizione iniziale fittizia)
% hold on;
% h_arm = plot3([0, 0], [0, 0], [0, 0], 'b', 'LineWidth', 4);   % Arto orizzontale
% h_pend = plot3([0, 0], [0, 0], [0, 0], 'r', 'LineWidth', 3);  % Asta del pendolo
% h_bob = plot3(0, 0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % Massa del pendolo
% hold off;
% 
% % Loop di animazione (aggiorna ogni 10 step per velocità)
% for k = 1:10:N
% 
%     theta1 = x_lin(1,k);  % Angolo dell'arto
%     theta2 = x_lin(2,k);  % Angolo del pendolo (misurato in rad)
% 
%     % Calcola la posizione dell'estremità dell'arto
%     arm_tip = [arm_length*cos(theta1), arm_length*sin(theta1), 0];
% 
%     % Calcola la posizione della massa del pendolo
%     phi =  pi - theta2; % Pendolo verticale verso l'alto quando theta2 = pi
%      pend_disp = [ -pendulum_length * cos(theta1) * sin(theta2), ...
%                    -pendulum_length * sin(theta1) * sin(theta2), ...
%                    -pendulum_length * cos(theta2) ];
%     bob_pos = arm_tip + pend_disp;
% 
%     % Aggiorna solo se gli oggetti sono validi
%     if isvalid(h_arm) && isvalid(h_pend) && isvalid(h_bob)
%         set(h_arm, 'XData', [0, arm_tip(1)], 'YData', [0, arm_tip(2)], 'ZData', [0, arm_tip(3)]);
%         set(h_pend, 'XData', [arm_tip(1), bob_pos(1)], 'YData', [arm_tip(2), bob_pos(2)], 'ZData', [arm_tip(3), bob_pos(3)]);
%         set(h_bob, 'XData', bob_pos(1), 'YData', bob_pos(2), 'ZData', bob_pos(3));
%     end
%     drawnow;
% end

% SIMULATION LQR ON THE REAL SYSTEM
x = zeros(5, N);
x(:,1) = x_lin(:,1);%[0; 5*pi/180; 0; 0; 0];  % Starting from the up position of the pendulum
    
% Optional disturbance vector
d = zeros(1, N);
d(1, N/2:N/2+1) = 0;

% Step-by-step simulation of the nonlinear system
y = x(:,1);  % Initialize current state
for k = 1:N-1
    V = max(min(-K*y,10),-10);                 % Control input (saturates at ±10)
    u2 = d(1, k);                  % Disturbance input (if any)
    
    % Simulate the nonlinear dynamics over one time step using ODE45
    [~, y_next] = ode45(@(t, y) furuta_nonlinear(y, V, u2), [0 Ts], y);
    
    % Update the state: take the final point from the ODE solution
    y = y_next(end, :)';
    x(:, k+1) = y;
end

% Plotting Results
figure(3);
subplot(2,2,1);
stairs(0:Ts:T_sim-Ts, x(1,:));
title('State x_1 (\theta_1)');
grid on;

subplot(2,2,2);
stairs(0:Ts:T_sim-Ts, x(2,:));
title('State x_2 (\theta_2)');
grid on;
yline(pi, '--r', 'Reference \pi');  % Reference line for theta2

subplot(2,2,3);
stairs(0:Ts:T_sim-Ts, x(3,:));
title('State x_3 (\theta_1 dot)');
grid on;

subplot(2,2,4);
stairs(0:Ts:T_sim-Ts, x(4,:));
title('State x_4 ({\theta}_2 dot)');
grid on;

xlabel('Time (s)');
sgtitle('State Evolution with LQR Control NOn Linear system');

%% Nonlinear Dynamics Function
function dydt = furuta_nonlinear(y, V, u2)
    % Unpack state variables
    theta1     = y(1);  % Angle of first link (rotor)
    theta2     = y(2);  % Angle of second link (pendulum)
    theta1_dot = y(3);  % Angular velocity of first link
    theta2_dot = y(4);  % Angular velocity of second link
    i_motor    = y(5);  % Motor current

    % System physical parameters
    L1 = 0.278;    % Length parameter 1 [m]
    L2 = 0.300;    % Length parameter 2 [m]
    l1 = 0.150;    % Center of mass distance for first link [m]
    l2 = 0.148;    % Center of mass distance for second link [m]
    m1 = 0.300;    % Mass of first link [kg]
    m2 = 0.075;    % Mass of second link [kg]
    J1 = 2.48e-2;  % Inertia of first link [kg*m^2]
    J2 = 3.86e-3;  % Inertia of second link [kg*m^2]
    b1 = 1.00e-4;  % Damping coefficient for first link
    b2 = 2.80e-4;  % Damping coefficient for second link
    g  = 9.81;     % Gravitational acceleration [m/s^2]
    Km = 0.090;    % Motor constant
    Lm = 0.005;    % Motor inductance [H]
    Rm = 7.80;     % Motor resistance [Ohm]
    
    % Effective inertia terms
    J0_hat = J1 + m1 * l1^2 + m2 * L1^2;
    J2_hat = J2 + m2 * l2^2;
    
    % Common denominator for the dynamic equations
    denom = J0_hat * J2_hat + J2_hat^2 * sin(theta2)^2 - m2^2 * L1^2 * l2^2 * cos(theta2)^2;
    u1 = Km*y(5);
    % Compute angular acceleration for theta1 (first link)
    theta1_ddot = ( - J2_hat * b1 * theta1_dot ...
                    - m2 * L1 * l2 * cos(theta2) * b2 * theta2_dot ...
                    - J2_hat^2 * sin(2 * theta2) * theta1_dot * theta2_dot ...
                    + 0.5 * J2_hat * m2 * L1 * l2 * cos(theta2) * sin(2 * theta2) * theta1_dot^2 ...
                    - J2_hat * m2 * L1 * l2 * sin(theta2) * theta2_dot^2 ...
                    + J2_hat * u1 ...
                    + m2 * L1 * l2 * cos(theta2) * u2 ...
                    + 0.5 * m2^2 * l2^2 * L1 * sin(2 * theta2) * g) / denom;
    
    % Compute angular acceleration for theta2 (second link)
    theta2_ddot = ( - m2 * L1 * l2 * cos(theta2) * b1 * theta1_dot ...
                    - b2 * (J0_hat + J2_hat * sin(theta2)^2) * theta2_dot ...
                    - m2 * L1 * l2 * J2_hat * cos(theta2) * sin(2 * theta2) * theta1_dot * theta2_dot ...
                    - 0.5 * sin(2 * theta2) * (J0_hat * J2_hat + J2_hat^2 * sin(theta2)^2) * theta1_dot^2 ...
                    - 0.5 * m2^2 * L1^2 * l2^2 * sin(2 * theta2) * theta2_dot^2 ...
                    + m2 * L1 * l2 * cos(theta2) * u1 ...
                    + (J0_hat + J2_hat * sin(theta2)^2) * u2 ...
                    + m2 * l2 * sin(theta2) * (J0_hat + J2_hat * sin(theta2)^2) * g) / denom;
    
    % Motor current dynamics
    i_motor_dot = (V-Rm*i_motor-Km*theta1_dot)/Lm;
    
    % Return state derivatives
    dydt = [theta1_dot; theta2_dot; theta1_ddot; theta2_ddot; i_motor_dot];
end
%% COMPARISON
% figure();
% subplot(2,2,1);
% stairs(0:Ts:T_sim-Ts, x(1,:), 'b', 'LineWidth', 1.5); hold on;
% stairs(0:Ts:T_sim-Ts, x_lin(1,:), 'r--', 'LineWidth', 1.5);
% title('State x_1 (\theta_1)');
% grid on;
% legend('Nonlinear', 'Linearized');
% hold off;
% 
% % Subplot for state x2 (theta2)
% subplot(2,2,2);
% stairs(0:Ts:T_sim-Ts, x(2,:), 'b', 'LineWidth', 1.5); hold on;
% stairs(0:Ts:T_sim-Ts, x_lin(2,:), 'r--', 'LineWidth', 1.5);
% title('State x_2 (\theta_2)');
% grid on;
% yline(pi, '--k', 'Reference \pi');
% legend('Nonlinear', 'Linearized');
% hold off;
% 
% % Subplot for state x3 (theta1 dot)
% subplot(2,2,3);
% stairs(0:Ts:T_sim-Ts, x(3,:), 'b', 'LineWidth', 1.5); hold on;
% stairs(0:Ts:T_sim-Ts, x_lin(3,:), 'r--', 'LineWidth', 1.5);
% title('State x_3 (\theta_1 dot)');
% grid on;
% legend('Nonlinear', 'Linearized');
% hold off;
% 
% % Subplot for state x4 (theta2 dot)
% subplot(2,2,4);
% stairs(0:Ts:T_sim-Ts, x(4,:), 'b', 'LineWidth', 1.5); hold on;
% stairs(0:Ts:T_sim-Ts, x_lin(4,:), 'r--', 'LineWidth', 1.5);
% title('State x_4 (\theta_2 dot)');
% grid on;
% legend('Nonlinear', 'Linearized');
% hold off;
% 
% xlabel('Time (s)');
% sgtitle('State Evolution Comparison: Nonlinear vs. Linearized Systems');
% 
% figure(5);
% stairs(0:Ts:T_sim-Ts, x(5,:), 'b', 'LineWidth', 1.5); hold on;
% stairs(0:Ts:T_sim-Ts, x_lin(5,:), 'r--', 'LineWidth', 1.5);
% title('State x_5 (Motor Current)');
% xlabel('Time (s)');
% ylabel('Current (A)');
% grid on;
% legend('Nonlinear', 'Linearized');
% hold off;
%% Graphic simulation
% % Creazione della figura
% fig = figure;
% axis equal; grid on;
% xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% view(3);
% xlim([-0.5, 0.5]); ylim([-0.5, 0.5]); zlim([-0.2, 0.5]);
% title('Simulazione Grafica del Furuta Pendulum con Controllore LQR');
% 
% % Parametri per la visualizzazione
% arm_length = L1;         % Lunghezza dell'arto orizzontale
% pendulum_length = L2;    % Lunghezza del pendolo
% 
% % Inizializza oggetti grafici (con posizione iniziale fittizia)
% hold on;
% h_arm = plot3([0, 0], [0, 0], [0, 0], 'b', 'LineWidth', 4);   % Arto orizzontale
% h_pend = plot3([0, 0], [0, 0], [0, 0], 'r', 'LineWidth', 3);  % Asta del pendolo
% h_bob = plot3(0, 0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % Massa del pendolo
% hold off;
% 
% % Loop di animazione (aggiorna ogni 10 step per velocità)
% for k = 1:10:N
% 
%     theta1 = x(1,k);  % Angolo dell'arto
%     theta2 = x(2,k);  % Angolo del pendolo (misurato in rad)
% 
%     % Calcola la posizione dell'estremità dell'arto
%     arm_tip = [arm_length*cos(theta1), arm_length*sin(theta1), 0];
% 
%     % Calcola la posizione della massa del pendolo
%     phi =  pi - theta2; % Pendolo verticale verso l'alto quando theta2 = pi
%     pend_disp = [ -pendulum_length * cos(theta1)*sin(phi), ...
%                   -pendulum_length * sin(theta1)*sin(phi), ...
%                    pendulum_length*cos(phi) ];
%     bob_pos = arm_tip + pend_disp;
% 
%     % Aggiorna solo se gli oggetti sono validi
%     if isvalid(h_arm) && isvalid(h_pend) && isvalid(h_bob)
%         set(h_arm, 'XData', [0, arm_tip(1)], 'YData', [0, arm_tip(2)], 'ZData', [0, arm_tip(3)]);
%         set(h_pend, 'XData', [arm_tip(1), bob_pos(1)], 'YData', [arm_tip(2), bob_pos(2)], 'ZData', [arm_tip(3), bob_pos(3)]);
%         set(h_bob, 'XData', bob_pos(1), 'YData', bob_pos(2), 'ZData', bob_pos(3));
%     end
%     pause(0);
%     drawnow;
% end
