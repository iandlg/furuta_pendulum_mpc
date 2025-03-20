clc; close all;
import mpt3.*

% Define system matrix (autonomous system)
A = LTI.Aclosed;  

% Define the autonomous system in MPT3
system = LTISystem('A', A);
system.x.min = [-1; -1; -1; -1; -1];
system.x.max = [1; 1; 1; 1; 1];

% Compute the Maximum Positively Invariant (MPI) set
%MPI_set = system.invariantSet();
MPI_set.isEmptySet()

% Plot the MPI set
figure;
plot(MPI_set);
title('Maximum Positively Invariant (MPI) Set');
xlabel('x_1'); ylabel('x_2');
grid on;

%% MPC implementation
clc; close all;
% 0. HyperParameters
A = [1 1 1 0; 0 1 0 1; 1 0 1 2; 1 1 1 1];  
%A = [1 1 1; 0 1 0; 1 0 1];  

B = [0; 1; 1; 1];
LTI.A = A_d; LTI.B = B_c;
Co = ctrb(LTI.A, LTI.B);
disp(['The rank of the controlability matrix of the linearized matrix pair (A_d B_d) is : ', num2str(rank(Co))])


options = sdpsettings('verbose',0,'solver','quadprog');
dim.N = 50;      % horizon
dim.nx = size(LTI.A,1);      % system order
dim.nu = 1;      % input order
param.T = 100;    % simulation number of steps
param.eps = 15*pi/180; % deviation from equilibrium

% 1. Define LQR weighting matrices
cost.Q = diag([0.000000001;0.00001;1;1;100]);  % Weighting on states (identity matrix)
cost.R = 100;             % Weighting on control input (scalar, as there's only one control input)

% Let Qf = P (solution to DARE)
[LTI.K, cost.Qf, ~] = dlqr(LTI.A, LTI.B, cost.Q, cost.R); % Optimal feedback gain

% Lifting the cost matrices for them to fit the optimization problem
cost.Ql = kron(eye(dim.N),cost.Q);cost.Ql = blkdiag(cost.Ql, cost.Qf);
cost.Rl = kron(eye(dim.N),cost.R);

% 2. Define lifted system dynamics
[LTI.Al, LTI.Bl] = predmodgen_state(LTI, dim);

% 3. Define constraints
%   State constraints : limit on theta 1 and 2 to stay linear
con.xmin = -5*ones(dim.nx,1);
con.xmax = 5*ones(dim.nx,1);
con.F = [eye(dim.nx); -eye(dim.nx)];
con.f = [con.xmax;  -con.xmin];

%   Input Constraint : set U st Gu \leq g
con.G = [1; -1]; con.g = [1; 1];

%   Terminal State Constraint : assume the same as state to be ok
LTI.K = -LTI.K;
LTI.Aclosed = LTI.A + LTI.B*LTI.K;

LTI.furutaClosed = LTISystem('A', LTI.Aclosed);
eig(LTI.Aclosed)

LTI.furutaClosed.x.min = con.xmin; 
LTI.furutaClosed.x.max = con.xmax;

% con.terminal_state_set = LTI.furutaClosed.invariantSet();
disp(['Terminal state set computed succesfully. Empty : ', num2str(con.terminal_state_set.isEmptySet())])
% tol1 = 0.2; tol2 = 0.2; tol3 = 1;
% con.Ff = con.F;
% con.ff = [tol1; tol1; tol2; tol2; tol3; tol3; tol3; tol3; tol3; tol3];

con.Ff = con.terminal_state_set.A; con.ff = con.terminal_state_set.b;
con.state_set = Polyhedron(con.F,con.f);
if con.state_set.contains(con.terminal_state_set)
    disp('P is fully contained in Q.');
else
    disp('P is NOT fully contained in Q.');
end

%   Lifted Constraints
con.Fl = kron(eye(dim.N), con.F); con.Fl = blkdiag(con.Fl, con.Ff);
con.fl = repmat(con.f,dim.N,1); con.fl = [con.fl; con.ff];

con.Gl = kron(eye(dim.N), con.G);
con.gl = repmat(con.g,dim.N,1);

% Initialize solver
u = sdpvar(repmat(dim.nu,1,dim.N),repmat(1,1,dim.N)); %#ok<*RPMT1>
x = sdpvar(repmat(dim.nx,1,dim.N+1),repmat(1,1,dim.N+1));

constraints = [];
objective = 0;
for k = 1:dim.N
 objective = objective + x{k}'*cost.Q*x{k} + u{k}'*cost.R*u{k};
 constraints = [constraints, x{k+1} == LTI.A*x{k} + LTI.B*u{k}];
 constraints = [constraints, -10 <= u{k}<= 10, -50<=x{k+1}<=50];
end
objective = objective + x{dim.N+1}'*cost.Qf*x{dim.N+1};

parameters_in = x{1};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','quadprog'),parameters_in,solutions_out);


% 4. Simulation : Receding horizon implementation for the constrained control problem
x = zeros(dim.nx,param.T);  % state vector with 0 at linearization point
x(:,1) = [0;param.eps;0;-1;0];  % starting from the up position of the pendulum with 

u_rec = zeros(dim.nu,param.T); % input vector
x_nonlin = zeros(dim.nx, param.T);
x_nonlin(:,1) = xref + x(:,1);
xref  = [0;pi;0;0;0];

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

    % Compute the state/output evolution             + LTI.Bdist*d(:,k)
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u_rec(:,k);

    % sim real siyst
    % Simulate the nonlinear dynamics over one time step using ODE45
    [~, y_next] = ode45(@(t, y) furuta_nonlinear(y, u_rec(:,k), 0), [0 Ts], x_nonlin(:,k));
    x_nonlin(:,k+1) = y_next(end,:)';
    clear u_con
end


% Plotting Results
figure(4);
subplot(3,2,1);
stairs(x(1,:));
title('State x_1 (\theta_1)');
grid on;

subplot(3,2,2);
stairs(x(2,:));
title('State x_2 (\theta_2)');
grid on;
yline(0, '--r', 'Reference \pi');  % Reference line for theta2

subplot(3,2,3);
stairs(x(3,:));
title('State x_3 (\theta_1 dot)');
grid on;

subplot(3,2,4);
stairs(x(4,:));
title('State x_4 ({\theta}_2 dot)');
grid on;

subplot(3,2,5);
stairs(x(5,:));
title('State x_5 (i)');
grid on;

subplot(3,3,6);
stairs(u_rec(1,:));
title('Input u (V)');
grid on;

xlabel('Time (s)');
sgtitle('State Evolution with MPC linear state control');

% Plotting Results
figure(5);
subplot(3,2,1);
stairs(x_nonlin(1,:));
title('State x_1 (\theta_1)');
grid on;

subplot(3,2,2);
stairs(x_nonlin(2,:));
title('State x_2 (\theta_2)');
grid on;
yline(0, '--r', 'Reference \pi');  % Reference line for theta2

subplot(3,2,3);
stairs(x_nonlin(3,:));
title('State x_3 (\theta_1 dot)');
grid on;

subplot(3,2,4);
stairs(x_nonlin(4,:));
title('State x_4 ({\theta}_2 dot)');
grid on;

subplot(3,2,5);
stairs(x_nonlin(5,:));
title('State x_5 (i)');
grid on;

subplot(3,3,6);
stairs(u_rec(1,:));
title('Input u (V)');
grid on;

xlabel('Time (s)');
sgtitle('State Evolution with MPC linear state control');

%%
yalmip('clear')
clear all

% Model data
A = [2 -1;1 0.2];
B = [1;0];
nx = 2; % Number of states
nu = 1; % Number of inputs

% MPC data
Q = eye(2);
R = 2;
N = 7;

% Initial state
x0 = [3;1];

u = sdpvar(repmat(nu,1,N),repmat(1,1,N)); %#ok<*RPMT1>
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

constraints = [];
objective = 0;
for k = 1:N
 objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
 constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
 constraints = [constraints, -1 <= u{k}<= 1, -5<=x{k+1}<=5];
end
objective = objective + x{N+1}'*Q*x{N+1};

parameters_in = x{1};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','quadprog'),parameters_in,solutions_out);

x = [3;1];
clf;
hold on
implementedU = [];
for i = 1:15
    inputs = {x};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end
  stairs(i:i+length(U)-1,U,'r')
  x = A*x + B*U(1);
  pause(0.05)
  stairs(i:i+length(U)-1,U,'k')
  implementedU = [implementedU;U(1)];
end
stairs(implementedU,'b')

%% 

yalmip('clear')
clear all

% Model data
A = [2 -1;1 0.2];
B = sdpvar(2,1);
E = [1;1];
nx = 2; % Number of states
nu = 1; % Number of inputs

% MPC data
Q = eye(2);
R = 2;
N = 20;


ny = 1;
C = [1 0];

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));
d = sdpvar(1);
pastu = sdpvar(1);

constraints = [-.1 <= diff([pastu u{:}]) <= .1];
objective = 0;
for k = 1:N
    objective = objective + (C*x{k}-r{k})'*(C*x{k}-r{k}) + u{k}'*u{k};
    constraints = [constraints, x{k+1} == A*x{k}+B*u{k}+E*d];
    constraints = [constraints, -1 <= u{k}<= 1, -5<=x{k+1}<=5];
end
objective = objective + (C*x{N+1}-r{N+1})'*(C*x{N+1}-r{N+1});

parameters_in = {x{1},[r{:}],d,pastu,B};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','quadprog'),parameters_in,solutions_out);
x = [0;0];
clf;
disturbance = randn(1)*.01;
oldu = 0;
hold on
xhist = x;
for i = 1:300
    if i < 50
        Bmodel = [1;0];
    else
        Bmodel = [.9;.1];
    end
    future_r = 3*sin((i:i+N)/40);    
    inputs = {x,future_r,disturbance,oldu,Bmodel};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};oldu = U(1);
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end    
    subplot(1,2,1);stairs(i:i+length(U)-1,U,'r')
    subplot(1,2,2);cla;stairs(i:i+N,X(1,:),'b');hold on;stairs(i:i+N,future_r(1,:),'k')
    stairs(1:i,xhist(1,:),'g')    
    x = A*x + Bmodel*U(1)+E*disturbance;
    xhist = [xhist x];
    pause(0.05)   
    % The measured disturbance actually isn't constant, it changes slowly
    disturbance = 0.99*disturbance + 0.01*randn(1);
end




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
    theta1_ddot = ( -J2_hat * b1 * theta1_dot ...
                    + m2 * L1 * l2 * cos(theta2) * b2 * theta2_dot ...
                    - J2_hat^2 * sin(2 * theta2) * theta1_dot * theta2_dot ...
                    - 0.5 * J2_hat * m2 * L1 * l2 * cos(theta2) * sin(2 * theta2) * theta1_dot^2 ...
                    + J2_hat * m2 * L1 * l2 * sin(theta2) * theta2_dot^2 ...
                    + J2_hat * u1 ...
                    - m2 * L1 * l2 * cos(theta2) * u2 ...
                    + 0.5 * m2^2 * l2^2 * L1 * sin(2 * theta2) * g) / denom;
    
    % Compute angular acceleration for theta2 (second link)
    theta2_ddot = ( m2 * L1 * l2 * cos(theta2) * b1 * theta1_dot ...
                    - b2 * (J0_hat + J2_hat * sin(theta2)^2) * theta2_dot ...
                    + m2 * L1 * l2 * J2_hat * cos(theta2) * sin(2 * theta2) * theta1_dot * theta2_dot ...
                    - 0.5 * sin(2 * theta2) * (J0_hat * J2_hat + J2_hat^2 * sin(theta2)^2) * theta1_dot^2 ...
                    - 0.5 * m2^2 * L1^2 * l2^2 * sin(2 * theta2) * theta2_dot^2 ...
                    - m2 * L1 * l2 * cos(theta2) * u1 ...
                    + (J0_hat + J2_hat * sin(theta2)^2) * u2 ...
                    - m2 * l2 * sin(theta2) * (J0_hat + J2_hat * sin(theta2)^2) * g) / denom;
    
    % Motor current dynamics
    i_motor_dot = (V-Rm*i_motor-Km*theta1_dot)/Lm;
    
    % Return state derivatives
    dydt = [theta1_dot; theta2_dot; theta1_ddot; theta2_ddot; i_motor_dot];
end