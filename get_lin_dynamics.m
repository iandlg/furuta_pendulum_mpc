function [Adisc,Bcontrol, Bdisturb] = get_lin_dynamics(y_eq, u_eq, Ts)
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

J0_hat = J1 + m1*l1^2 + m2*L1^2;
J2_hat = J2 + m2*l2^2;

% Define symbolic variables for linearization
syms y1 y2 y3 y4 y5 V u2 real
dyn.y = [y1; y2; y3; y4; y5];  % State vector [theta1, theta2, theta1_dot, theta2_dot]
dyn.u = [V; u2];  % Input vector [u1, u2]

% Equations of Motion (Nonlinear)
    dy1 = y3;
    dy2 = y4;
    dy3 = (- J2_hat*b1*y3 ...
           - m2*L1*l2*cos(y2)*b2*y4 ...
           - J2_hat^2*sin(2*y2)*y3*y4 ...
           + 0.5*J2_hat*m2*L1*l2*cos(y2)*sin(2*y2)*y3^2 ...
           - J2_hat*m2*L1*l2*sin(y2)*y4^2 ...
           + J2_hat*Km*y5 ...
           + m2*L1*l2*cos(y2)*u2 ...
           + 0.5*m2^2*l2^2*L1*sin(2*y2)*g)... 
           /(J0_hat*J2_hat + J2_hat^2*sin(y2)^2 - m2^2*L1^2*l2^2*cos(y2)^2);
    dy4 = ( - m2*L1*l2*cos(y2)*b1*y3 ...
            - b2*(J0_hat + J2_hat*sin(y2)^2)*y4 ...
            - m2*L1*l2*J2_hat*cos(y2)*sin(2*y2)*y3*y4 ...
            - 0.5*sin(2*y2)*(J0_hat*J2_hat + J2_hat^2*sin(y2)^2)*y3^2 ...
            - 0.5*m2^2*L1^2*l2^2*sin(2*y2)*y4^2 ...
            + m2*L1*l2*cos(y2)*Km*y5 ...
            + (J0_hat + J2_hat*sin(y2)^2)*u2...
            + m2*l2*sin(y2)*(J0_hat + J2_hat*sin(y2)^2)*g)...
            / (J0_hat*J2_hat + J2_hat^2*sin(y2)^2 - m2^2*L1^2*l2^2*cos(y2)^2);
    dy5 = (V - Rm*y5 - Km*y3)/Lm;
% System of equations
dyn.f = [dy1; dy2; dy3; dy4; dy5];


% Linearize the system around equilibrium point
% Jacobian of the system with respect to state variables and inputs
A = jacobian(dyn.f, dyn.y);  % Jacobian with respect to the state vector y
B = jacobian(dyn.f, dyn.u);  % Jacobian with respect to the input vector u

% Evaluate the Jacobians at the equilibrium point
A_eq = double(subs(A, [dyn.y; dyn.u], [y_eq; u_eq]));  % Evaluate A at equilibrium
B_eq = double(subs(B, [dyn.y; dyn.u], [y_eq; u_eq]));  % Evaluate B at equilibrium

% Display the linearized system matrices
% disp('Linearized A Matrix:');
% disp(A_eq);
% disp('Linearized B Matrix:');
% disp(B_eq);

% DISCRETIZATION
% Discretization using zero-order hold (ZOH)
sys_disc = c2d(ss(A_eq, B_eq, eye(5), zeros(5, 2)), Ts, 'zoh');

% Extract discrete-time system matrices
Adisc = sys_disc.A;
Bdisc = sys_disc.B;

Bcontrol = Bdisc(:,1);
Bdisturb = Bdisc(:,2);

end