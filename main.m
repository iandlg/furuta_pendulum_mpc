% Futura Pendulum Simulation (Nonlinear Dynamics with Inputs)

clear; clc; close all;
eps = 0.08; % deviation from initial condition
dist = +0.1; %intensity of the disturbance torque
%% Parameters (Original System)
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

%% Define Control Inputs (Torque Inputs)
V = @(t) -10;  % Voltage applied to the motor (V)
u2 = @(t) 0;  % External disturbance torque on Arm 2 (Nm)

%% Equations of Motion with Inputs
 pendulum_dynamics = @(t, y) [
        y(3);
        y(4);
        (- J2_hat*b1*y(3) ...
         - m2*L1*l2*cos(y(2))*b2*y(4) ...
         - J2_hat^2*sin(2*(y(2)))*y(3)*y(4) ...
         + 0.5*J2_hat*m2*L1*l2*cos(y(2))*sin(2*(y(2)))*y(3)^2 ...
         - J2_hat*m2*L1*l2*sin((y(2)))*y(4)^2 ...
         + J2_hat*Km*y(5) ...
         + m2*L1*l2*cos((y(2)))*u2(t) ...
         + 0.5*m2^2*l2^2*L1*sin(2*(y(2)))*g) ...
         /(J0_hat*J2_hat + J2_hat^2*sin(y(2))^2-m2^2*L1^2*l2^2*cos(y(2))^2);
    
        (- m2*L1*l2*cos(y(2))*b1*y(3) ...
         - b2*(J0_hat+J2_hat*sin(y(2))^2)*y(4) ...
         - m2*L1*l2*J2_hat*cos(y(2))*sin(2*(y(2)))*y(3)*y(4) ...
         - 0.5*sin(2*(y(2)))*(J0_hat*J2_hat + J2_hat^2*sin(y(2))^2)*y(3)^2 ...
         - 0.5*m2^2*L1^2*l2^2*sin(2*(y(2)))*y(4)^2 ...
         + m2*L1*l2*cos(y(2))*Km*y(5) ...
         + (J0_hat+J2_hat*sin(y(2))^2)*u2(t) ...
         + m2*l2*sin(y(2))*(J0_hat+J2_hat*sin(y(2))^2)*g)...
         /(J0_hat*J2_hat + J2_hat^2*sin(y(2))^2-m2^2*L1^2*l2^2*cos(y(2))^2);
    
         (V(t)-Rm*y(5)-Km*y(3))/Lm;
    ];
%% Simulation Settings
Tspan = [0 4];
Y0 = [0; 5*pi/180; 0; 0; 0]; % Initial conditions: [theta1, theta2, theta1_dot, theta2_dot]

%% Solve ODE using a more stable solver
[t, Y] = ode15s(pendulum_dynamics, Tspan, Y0);

%% Plot Results
figure(1);
plot(t, Y(:,1), 'r', t, Y(:,2), 'b');
legend('\theta_1 (Arm Rotation)', '\theta_2 (Pendulum Angle)');
xlabel('Time (s)'); ylabel('Angle (rad)'); title('Futura Pendulum Simulation with Input Torque');

grid on;
%% LINEARIZATION OF THE SYSTEM 
%% Define symbolic variables for linearization
syms y1 y2 y3 y4 y5 V u2 real
dyn.y = [y1; y2; y3; y4; y5];  % State vector [theta1, theta2, theta1_dot, theta2_dot]
dyn.u = [V; u2];  % Input vector [u1, u2]

%% Equations of Motion (Nonlinear)
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


%% MPC implementation - Linearized
clc; close all;
disp("MPC implementation - Linearized")
% Get linearized system
param.y_eq = [0; 0; 0; 0; 0];  % Equilibrium point [theta1, theta2, theta1_dot, theta2_dot, i_motor]
param.u_eq = [0; 0];  % Equilibrium input [u1, u2]
param.Ts = 0.2;

[LTI.A, LTI.B, LTI.Bdist] = get_lin_dynamics(dyn,param.y_eq,param.u_eq,param.Ts);

% Parameters
options = sdpsettings('verbose',0,'solver','quadprog');
dim.N =5;      % horizon
dim.nx = size(LTI.A,1);      % system order
dim.nu = 1;      % input order
dim.nd = 1;  %number of disturbance
param.time = 0:param.Ts:30;
param.T = length(param.time);    % simulation number of steps
param.eps = 4.8*pi/180; % deviation from equilibrium

Co = ctrb(LTI.A, LTI.B);
disp(['The rank of the controlability matrix of the linearized matrix pair (A_d B_d) is : ', num2str(rank(Co))])

% 1. Define LQR weighting matrices
cost.Q = diag([10;10;1;1;1]);  % Weighting on states (identity matrix)
cost.R = 15;             % Weighting on control input (scalar, as there's only one control input)

% Let Qf = P (solution to DARE)
[LTI.K, cost.Qf, ~] = dlqr(LTI.A, LTI.B, cost.Q, cost.R); % Optimal feedback gain

% 3. Define constraints
%   State constraints : limit on theta 1 and 2 to stay linear
con.xmax = [pi; 4.8*pi/180; inf; inf; inf];
con.xmin = -con.xmax;

%   Input Constraint : set U st Gu \leq g
con.umax = 10;
con.umin = -con.umax;
LTI.x0 = [0;param.eps;0;0;0];
%   Terminal State Constraint : assume the same as state to be ok
[con.Ff, con.ff] = get_term_state_constraints(LTI, con, dim);

% Initialize solver
u = sdpvar(repmat(dim.nu,1,dim.N),ones(1,dim.N)); 
x = sdpvar(repmat(dim.nx,1,dim.N+1),ones(1,dim.N+1));

constraints = [];
objective = 0;
x_ref = [0;0;0;0;0];
for k = 1:dim.N
 objective = objective + (x{k}-x_ref)'*cost.Q*(x{k}-x_ref) + u{k}'*cost.R*u{k};
 constraints = [constraints, x{k+1} == LTI.A*x{k} + LTI.B*u{k}];
 constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmin <= x{k+1}<= con.xmax];
end
% constraints = [constraints; con.Ff*x{dim.N+1} <= con.ff];
% constraints = [constraints; con.xmin <= x{dim.N+1}<= con.xmax];

objective = objective + 2*(x{dim.N+1}-x_ref)'*cost.Qf*(x{dim.N+1}-x_ref);

parameters_in = x{1};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','quadprog'),parameters_in,solutions_out);


% 4. Simulation : Receding horizon implementation for the constrained control problem
x = zeros(dim.nx,param.T);  % state vector with 0 at linearization point
x(:,1) = [0;param.eps;0;0;0];  % starting from the up position of the pendulum with 
d = zeros(1, param.T); d(:,floor(param.T/2)) = 0.01;
u_rec = zeros(dim.nu,param.T); % input vector

for k=1:param.T-1
    inputs = {x(:,k)};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    % Select the first input only
    u_rec(:,k) = U(1);
    % u_rec(:,k) = -LTI.K*x(:,k);
    % Compute the state/output evolution
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u_rec(:,k); % + LTI.Bdist*d(:,k)
end


plot_state(x, u_rec, 'State and Input Evolution with MPC linear state control');
% pause(10);clc;close all;
%% MPC implementation - Non Linear dynamics; full state knowledge
clc;
disp("MPC implementation - Non Linear dynamics; full state knowledge")
u_rec = zeros(dim.nu,param.T); % input vector
x_nonlin = zeros(dim.nx, param.T);
xref  = [0;0;0;0;0];

x_nonlin(:,1) = xref + [0;param.eps;0;0;0];

for k=1:param.T
    inputs = {x_nonlin(:,k)-xref};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    % Select the first input only
    u_rec(:,k) = U(1);
    % u_rec(:,k) = -LTI.K*x(:,k);


    % sim real siyst
    % Simulate the nonlinear dynamics over one time step using ODE45
    [~, y_next] = ode45(@(t, y) futura_nonlinear(y, u_rec(:,k), 0), [0 param.Ts], x_nonlin(:,k));
    x_nonlin(:,k+1) = y_next(end,:)';
end

plot_state(x_nonlin, u_rec, 'State and Input Evolution with MPC nonlinear state control');
pause(10);clc;close all;
%% MPC with observer, partial state knoledge
close all;clc;
disp("MPC implementation - Linear dynamics; Observer")

param.eps = 2*pi/180;
yref = [0;0];
x0 = [0;param.eps;0;0;0];
LTI.x0 = x0;

dim.nd = 1;
LTI.C = [1 0 0 0 0;
         0 1 0 0 0]; % assume access to theta1 and 2 angles

dim.ny = size(LTI.C, 1);
LTI.Cdist = zeros(dim.ny,dim.nd);
%{
A_d = 0.000;      % Ampiezza del disturbo
omega_d = 1 * pi * 0.5;  % Frequenza in rad/s

S = [0 omega_d; -omega_d 0];

% Matrice di transizione discreta
Sin_disc = expm(S * param.Ts);

LTI.Bdist = [LTI.Bdist(:,1), zeros(dim.nx,1)];
LTIe.A = [LTI.A, LTI.Bdist; 
          zeros(dim.nd, dim.nx), Sin_disc]; % Ampiezza inclusa
%}
LTIe.A = [LTI.A, LTI.Bdist; 
          zeros(dim.nd, dim.nx), eye(dim.nd)]; % Ampiezza inclusa
LTIe.B = [LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C = [LTI.C, LTI.Cdist];
% LTIe.x0 = [x0; A_d; 0];
LTIe.x0 = [x0; 0];

LTIe.yref = yref;

dime.nx=dim.nx+dim.nd;     %state dimension
dime.nu=1;     %input dimension
dime.ny=2;     %output dimension
dime.N=5;      %horizon

weighte.Q=blkdiag(diag([10,10,10,10,10]),zeros(dim.nd));            %weight on output
weighte.R=cost.R;                                   %weight on input
weighte.P=blkdiag(cost.Qf,zeros(dim.nd));  

xe=zeros(dime.nx,param.T+1);
y=zeros(dime.ny,param.T+1);
u_rec=zeros(dime.nu,param.T);
xehat=zeros(dime.nx,param.T+1);

xe(:,1)=LTIe.x0;
xehat(:,1)=[0; 0; 0; 0; 0; 0];

measurement_noise = 0.00*normrnd(0,1,dim.ny,1);
y(:,1)=LTIe.C*LTIe.x0 +  measurement_noise;


L = place(LTIe.A',LTIe.C',[0.6; 0.55; 0.5;0.65 ;0.7; 0.4])';

options = sdpsettings('verbose',0,'solver','quadprog');

con.xmax = [inf; inf*pi/180; inf; inf; inf];
con.xmin = -con.xmax;
con.xmaxe = [inf; 5*pi/180; inf; inf; inf; 10];
con.xmine = -con.xmaxe;
con.umax = 10;
con.umin = -con.umax;

[con.Ff, con.ff] = get_term_state_constraints(LTI, con, dim);

% optimizer for the controller
u = sdpvar(repmat(dime.nu,1,dime.N),ones(1,dime.N)); 
x = sdpvar(repmat(dime.nx,1,dime.N+1),ones(1,dime.N+1));
xr = sdpvar(dime.nx,1);
ur = sdpvar(dime.nu, 1);

constraints = [];
objective = 0;
for k = 1:dime.N
 objective = objective + (x{k}-xr)'*weighte.Q*(x{k}-xr) + (u{k}-ur)'*weighte.R*(u{k}-ur);
 constraints = [constraints, x{k+1} == LTIe.A*x{k} + LTIe.B*u{k}];
 constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmine <= x{k+1}<= con.xmaxe];
end
constraints = [constraints; con.Ff*x{dim.N+1}(1:5) <= con.ff];

objective = objective + (x{dime.N+1}-xr)'*weighte.P*(x{dime.N+1}-xr);

parameters_in = {x{1}, xr, ur};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,options,parameters_in,solutions_out);
%
%optimizer for the target selection
%{
xref = sdpvar(dim.nx,1);
uref = sdpvar(dim.nu,1);
dtilde = sdpvar(dim.nd,1);
constraints = [];
constraints = [
                [eye(dim.nx)-LTI.A -LTI.B;
                LTI.C zeros(dim.ny,dim.nu)]*[xref;uref]==[LTI.Bdist*dtilde; yref-LTI.Cdist*dtilde];
                con.xmin <= xref <= con.xmax;
                con.umin <= uref <= con.umax
               ];
objective = [];
% objective = uref'*uref;
objective = 0;
solutions_input = {dtilde};
solutions_output = {xref, uref};

target_selector = optimizer(constraints, objective, options, solutions_input, solutions_output);
%}
param.T = 150;
for k=1:param.T
    k
    if(k>=90&& k<=95 )
        xe(end,k) = 0.0002;
    else
        xe(end,k) = 0;
    end
    xe_0=xe(:,k);  
    dhat=xehat(end-dim.nd+1:end,k);
    
    %{
    % Compute optimal ss (online, at every iteration)
    inputs = {dhat};
    [solutions, diagnostics] = target_selector{inputs};
    xref = solutions{1};
    uref = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible for target selector');
    end
    xre = [xref; dhat];
    %}
    
    [ineqconstraint, eqconstraints]= constraintsgen(LTI, dim, dhat, yref,con);
    [xref, uref] = optimalss(LTI, dim, cost, ineqconstraint, eqconstraints);
    xre = [xref; dhat];
    
     inputs = {xehat(:,k), xre, uref};
     [solutions,diagnostics] = controller{inputs};    
     U = solutions{1};
     X = solutions{2};
     % if diagnostics == 1
     %     error('The problem is infeasible for controller');
     % end     

    % Select the first input only
    % u_rec(:,k)=U(1:dim.nu);
    u_rec(:,k)= -LTI.K*xehat(1:dim.nx,k);
    % Compute the state/output evolution
    process_noise = 0 * normrnd(0,1, [dime.nx, 1]);
    xe(:,k+1)=LTIe.A*xe_0 + LTIe.B*u_rec(:,k) + process_noise;
    measurement_noise = 0*normrnd(0,1,dim.ny,1);
    y(:,k+1)=LTIe.C*xe(:,k+1)+ measurement_noise;
        
    % Update extended-state estimation

    xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));

end
%
plot_state(xe , u_rec , 'State Evolution with output MPC observer',xehat);

%% Output MPC - with input and measurement noise, Kalman fiklter
disp("Output MPC - with input and measurement noise, Kalman filter")

param.eps = 1*pi/180; % deviation from equilibrium

Ob = obsv(LTI.A, LTI.C);
disp(['The rank of the observability matrix of the linearized matrix pair (A C) is : ', num2str(rank(Ob))])

% Kalman Filter variables
cov.pos = eye(dime.nx);
cov.measurement = eye(dime.ny);
cov.motion = 0.1*eye(dime.nx);

% Initialize solver

% Simulation : Receding horizon implementation for the constrained control problem
xhat = zeros(dime.nx,param.T);  % state vector with 0 at linearization point
xhat(:,1) = [0;0;0;0;0;0];  % starting from the up position of the pendulum with
x = zeros(dime.nx,param.T);
x(:,1) = [0;param.eps;0;0;0;0];       % real starting position for state vector
u_rec = zeros(dim.nu,param.T); % input vector

for k=1:param.T
    % Do measurement update 
    measurement_noise = 0.001*normrnd(0,1,dime.ny,1);
    y = LTIe.C*x(:,k) + measurement_noise;
    [xhat(:,k), cov.pos] = measurement_update(xhat(:,k), y, LTIe, cov);
    
    if(k>=90&& k<=95 )
        x(end,k) = 0.00;
    else
        x(end,k) = 0;
    end

    dhat=xehat(end-dim.nd+1:end,k);

    [ineqconstraint, eqconstraints]= constraintsgen(LTI, dim, dhat, yref,con);
    [xref, uref] = optimalss(LTI, dim, cost, ineqconstraint, eqconstraints);
    xre = [xref; dhat];

    % Compute optimal control action
    inputs = {xhat(:,k), xre, uref};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible for controller');
    end   
    
    % Select the first input only
    u_rec(:,k) = U(1);

    % Compute the state/output evolution
    process_noise = zeros(dime.nx, 1); process_noise(end,1) = 0.5* normrnd(0,1);
    x(:,k+1) = LTIe.A*x(:,k) + LTIe.B*u_rec(:,k) + process_noise;
    
    % Do dynamic update for next k
    [xhat(:,k+1), cov.pos] = dynamic_update(xhat(:,k), u_rec(:,k), LTIe, cov);
    
end
% Note we apply uniform noise to the measurement and noise only to the
% input


plot_state(x, u_rec, 'State Evolution with output MPC on linearized system,Kalman filter',xhat);

%% state MPC - reference tracking (partial rotation of base arm)
clc; close all;
disp("state MPC - reference tracking (partial rotation of base arm)")
% Get linearized system
param.y_eq = [0; 0; 0; 0; 0];  % Equilibrium point [theta1, theta2, theta1_dot, theta2_dot, i_motor]
param.u_eq = [0; 0];  % Equilibrium input [u1, u2]
param.Ts = 0.1;


[LTI.A, LTI.B, LTI.Bdist] = get_lin_dynamics(dyn,param.y_eq,param.u_eq,param.Ts);
% Parameters
dim.N = 5;      % horizon
dim.nx = size(LTI.A,1);      % system order
dim.nu = 1;      % input order
param.time = 0:param.Ts:50;     % simulation time
param.T = length(param.time);    % simulation number of steps

param.eps = 4*pi/180; % deviation from equilibrium

Co = ctrb(LTI.A, LTI.B);
disp(['The rank of the controlability matrix of the linearized matrix pair (A_d B_d) is : ', num2str(rank(Co))])

% 1. Define LQR weighting matrices
cost.Q = diag([1;1;1;1;1]);  % Weighting on states (identity matrix)
cost.R = 1;             % Weighting on control input (scalar, as there's only one control input)

% Let Qf = P (solution to DARE)
[LTI.K, cost.Qf, ~] = dlqr(LTI.A, LTI.B, cost.Q, cost.R); % Optimal feedback gain

% Initialize solver
options = sdpsettings('verbose',0,'solver','quadprog');

con.xmax = [10; 10*pi/180; 10; 10; 100];
con.xmin = -con.xmax;
con.umax = 100;
con.umin = -con.umax;

[con.Ff, con.ff] = get_term_state_constraints(LTI, con, dim);

u = sdpvar(repmat(dim.nu,1,dim.N),ones(1,dim.N)); 
x = sdpvar(repmat(dim.nx,1,dim.N+1),ones(1,dim.N+1));
xr = sdpvar(dim.nx,1);
ur = sdpvar(dim.nu, 1);

constraints = [];
objective = 0;
for k = 1:dim.N
 objective = objective + (x{k}-xr)'*cost.Q*(x{k}-xr) + (u{k}-ur)'*cost.R*(u{k}-ur);
 constraints = [constraints, x{k+1} == LTI.A*x{k} + LTI.B*u{k}];
 constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmin <= x{k+1}<= con.xmax];
end
constraints = [constraints; con.Ff*x{dim.N+1} <= con.ff];
% constraints = [constraints; con.xmin <= x{dim.N+1}<= con.xmax];

objective = objective + (x{dim.N+1}-xr)'*cost.Qf*(x{dim.N+1}-xr);

parameters_in = {x{1}, xr, ur};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,options,parameters_in,solutions_out);
%{
% Target selector
options = sdpsettings('verbose',0,'solver','quadprog');

xr = sdpvar(dim.nx,1);
ur = sdpvar(dim.nu,1);
yr = sdpvar(dim.ny,1);

constraints = [
                [eye(dim.nx)-LTI.A, -LTI.B;
                LTI.C, zeros(dim.ny, dim.nu)]*[xr;ur]==[zeros(dim.nx,1); yr];
                con.xmin <= xr <= con.xmax;
                con.umin <= ur <= con.umax
               ];
objective = ur'*ur;     % squared norm of u_ref to draw it down to 0

parameters_in = {yr};
solutions_out = {xr, ur};

target_selector = optimizer(constraints, objective, options, parameters_in, solutions_out);
%}
% Generate y_ref
y_ref = generate_reference(param.time, 360*pi/180,false) - param.y_eq(1:2);

% Simulation : Receding horizon implementation for the constrained control problem
xhat = zeros(dim.nx,param.T);  % state vector with 0 at linearization point
xhat(:,1) = [0;0;0;0;0];  % starting from the up position of the pendulum with
x = zeros(dim.nx,param.T);
x(:,1) = [0;param.eps;0;0;0];       % real starting position for state vector
u_rec = zeros(dim.nu,param.T); % input vector

for k=1:param.T-1
    % Do measurement update 
    measurement_noise = 0.001*normrnd(0,1,2,1);
    y = LTI.C*x(:,k)+measurement_noise;
    [xhat(:,k), cov.pos] = measurement_update(xhat(:,k), y, LTI, cov);
    %{
    % Target selection
    inputs = {y_ref(:,k)};
    [solutions, diagnostics] = target_selector{inputs};
    x_ref = solutions{1};
    u_ref = solutions{2};
    if diagnostics == 1
        error(['The target selection is infeasible (k = ',num2str(k),')' ]');
    end
    %}

    [ineqconstraint, eqconstraints] = constraintsgen_reference(LTI, dim, 0, y_ref(:,k),con);
    [x_ref, u_ref] = optimalss(LTI, dim, cost, ineqconstraint, eqconstraints);
    
    % Solve problem
    inputs = {xhat(:,k), x_ref, u_ref};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error(['The problem is infeasible (k = ',num2str(k),')' ]');
    end
    
    % Select the first input only
    u_rec(:,k) = U(1);

    % Compute the state/output evolution
    process_noise = zeros(dim.nx, 1); process_noise(end,1) = 0.5* normrnd(0,1);
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u_rec(:,k) + process_noise; % + LTI.Bdist*d(:,k)

    % Do dynamic update for next k
    [xhat(:,k+1), cov.pos] = dynamic_update(xhat(:,k), u_rec(:,k), LTI, cov);
end


% Plotting Results
figure(6); clf;
subplot(3,2,1);

plot(param.time,y_ref(1,:)); hold on;
stairs(param.time,xhat(1,:));

plot(param.time,x(1,:),'--', color='black');

title('State x_1 (\theta_1)');
grid on;

subplot(3,2,2);
plot(param.time,x(2,:),'--', color='black'); hold on;
stairs(param.time,xhat(2,:));
title('State x_2 (\theta_2)');
grid on;
yline(0, '--r', 'Reference \pi');  % Reference line for theta2

subplot(3,2,3);
plot(param.time,x(3,:),'--', color='black'); hold on;
stairs(param.time,xhat(3,:));
title('State x_3 (\theta_1 dot)');
grid on;

subplot(3,2,4);
plot(param.time,x(4,:),'--', color='black'); hold on;
stairs(param.time,xhat(4,:));
title('State x_4 ({\theta}_2 dot)');
grid on;

subplot(3,2,5);
plot(param.time,x(5,:),'--', color='black'); hold on;
stairs(param.time,xhat(5,:));
title('State x_5 (i)');
grid on;

subplot(3,2,6);
stairs(param.time,u_rec(1,:)); 
title('Input u (V)');
grid on;

xlabel('Time (s)');
sgtitle('State Evolution with MPC linear state control');
%% Adaptive MPC - state knowledge - no terminal set
clc; close all;

% Parameters
options = sdpsettings('verbose',0,'solver','quadprog');
dim.N = 50;      % horizon
dim.nx = size(LTI.A,1);      % system order
dim.nu = 1;      % input order
param.T = 100;    % simulation number of steps
param.eps = 1*pi/180; % deviation from equilibrium

Co = ctrb(LTI.A, LTI.B);
disp(['The rank of the controlability matrix of the linearized matrix pair (A_d B_d) is : ', num2str(rank(Co))])

% 1. Define LQR weighting matrices
cost.Q = diag([0.1;0.1;1;1;10]);  % Weighting on states (identity matrix)
cost.R = 1;             % Weighting on control input (scalar, as there's only one control input)

% Let Qf = P (solution to DARE)
[LTI.K, cost.Qf, ~] = dlqr(LTI.A, LTI.B, cost.Q, cost.R); % Optimal feedback gain

% Initialize solver
u = sdpvar(repmat(dim.nu,1,dim.N),ones(1,dim.N)); 
x = sdpvar(repmat(dim.nx,1,dim.N+1),ones(1,dim.N+1));
A = sdpvar(dim.nx, dim.nx);
B = sdpvar(dim.nx, dim.nu);

constraints = [];
objective = 0;
for k = 1:dim.N
 objective = objective + x{k}'*cost.Q*x{k} + u{k}'*cost.R*u{k};
 constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
 constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmin <= x{k+1}<= con.xmax];
end
% constraints = [constraints; con.Ff*x{dim.N+1} <= con.ff];
constraints = [constraints; con.xmin <= x{dim.N+1}<= con.xmax];

objective = objective + x{dim.N+1}'*cost.Qf*x{dim.N+1};

parameters_in = {x{1}, A, B};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','quadprog'),parameters_in,solutions_out);

% 4. Simulation : Receding horizon implementation for the constrained control problem
x = zeros(dim.nx,param.T);  % state vector with 0 at linearization point
x(:,1) = [0;param.eps;0;0;0];  % starting from the up position of the pendulum with 
u_rec = zeros(dim.nu,param.T); % input vector

for k=1:param.T
    % System linearization 
    [LTI.A, LTI.B, LTI.Bdisc] = get_lin_dynamics(dyn,x(:,k),[0;0],param.Ts);

    % Solve problem
    inputs = {x(:,k),LTI.A, LTI.B};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    % Select the first input only
    u_rec(:,k) = U(1);

    % Compute the state/output evolution
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u_rec(:,k); % + LTI.Bdist*d(:,k)
end
