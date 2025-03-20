function [Adisc,Bcontrol, Bdisturb] = get_lin_dynamics(dyn, y_eq, u_eq, Ts)
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
