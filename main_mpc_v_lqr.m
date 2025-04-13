%% MPC implementation - linearized - MPC v LQR
clc;
disp("MPC implementation - linearized - MPC v LQR")

% Get linearized system
param.y_eq = [0; 0; 0; 0; 0];  % Equilibrium point [theta1, theta2, theta1_dot, theta2_dot, i_motor]
param.u_eq = [0; 0];  % Equilibrium input [u1, u2]
param.Ts = 0.2;

[LTI.A, LTI.B, LTI.Bdist] = get_lin_dynamics(param.y_eq,param.u_eq,param.Ts);

% Simulation parameters
options = sdpsettings('verbose',0,'solver','quadprog');
dim.N = 10;      % horizon
dim.nx = size(LTI.A,1);      % system order
dim.nu = 1;     % input order
dim.nd = 1;     % number of disturbance
param.time = 0:param.Ts:10;
param.T = length(param.time);    % simulation number of steps
param.eps = 1*pi/180; % deviation from equilibrium -> max is 4.968

Co = ctrb(LTI.A, LTI.B);
disp(['The rank of the controlability matrix of the linearized matrix pair (A_d B_d) is : ', num2str(rank(Co))])

% 1. Define LQR weighting matrices
cost.Q = diag([0.01;100;0.1;10;1]);  % Weighting on states (identity matrix)
cost.R = 1;             % Weighting on control input (scalar, as there's only one control input)

% Let Qf = P (solution to DARE)
[LTI.K, cost.Qf, ~] = dlqr(LTI.A, LTI.B, cost.Q, cost.R); % Optimal feedback gain
cost.Qf = 2*cost.Qf;
cost.maxP = eye(size(cost.Qf,1))*max(eig(cost.Qf));
cost.maxP_10 = 10*eye(size(cost.Qf,1))*max(eig(cost.Qf));
cost.minQ = eye(size(cost.Q,1))*min(eig(cost.Q));

Vf = @(x)  1/2*x' * cost.Qf * x;

ell = @(x,u) 1/2*x' * cost.Q * x + 1/2*u' * cost.R * u;

% 3. Define constraints
%   State constraints : limit on theta 1 and 2 to stay linear
con.xmax = [1e2; 10*pi/180; 1e2; 1e2; 1e2];
con.xmin = -con.xmax;

%   Input Constraint : set U st Gu \leq g
con.umax = 10;
con.umin = -con.umax;
LTI.x0 = [0;param.eps;0;0;0];

load("mpt3_controller_data.mat")

%   Terminal State Constraint : assume the same as state to be ok
con.Ff = Tset_mpt3.A; con.ff = Tset_mpt3.b;


u = sdpvar(repmat(dim.nu,1,dim.N),ones(1,dim.N)); 
x = sdpvar(repmat(dim.nx,1,dim.N+1),ones(1,dim.N+1));

constraints = [];
objective = 0;
x_ref = [0;0;0;0;0];
for k = 1:dim.N
    objective = objective + 1/2*(x{k}-x_ref)'*cost.Q*(x{k}-x_ref) + 1/2*u{k}'*cost.R*u{k};
    constraints = [constraints, x{k+1} == LTI.A*x{k} + LTI.B*u{k}];
    constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmin <= x{k+1}<= con.xmax];
end
constraints = [constraints; con.Ff*x{dim.N+1} <= con.ff];
% constraints = [constraints; con.xmin <= x{dim.N+1}<= con.xmax];

objective = objective + 1/2*(x{dim.N+1}-x_ref)'*cost.Qf*(x{dim.N+1}-x_ref);

parameters_in = x{1};
solutions_out = {[u{:}], [x{:}],objective};

controller = optimizer(constraints, objective,options,parameters_in,solutions_out);


% 4. Simulation : Receding horizon implementation for the constrained control problem
x = zeros(dim.nx,param.T);  % state vector with 0 at linearization point
x(:,1) =[0;param.eps;0;0;0];  % starting from the up position of the pendulum with
% param.eps = -0.3316;
% [0;-0.1745;0;1.07;0];
% [0;-0.3316;0;1.9;0];
% [0;param.eps;0;0;0];

%x(:,1) = LTI.x0;
x_lqr = x;
d = zeros(1, param.T); d(:,floor(param.T/2)) = 0.01;
u_rec = zeros(dim.nu,param.T); % input vector
u_lqr = u_rec;
mpc_in_tset = zeros(1,param.T);
lqr_in_tset = mpc_in_tset;
a = zeros(1,param.T);       % Vf(x+) - Vf(x)
b = zeros(1,param.T);       % - l(x,u)

% cost function to verify the asusmptions a posteriori
obj = zeros(size(param.T));
alpha_f = zeros(size(param.T));
alpha_f_10 = zeros(size(param.T));
V_f = zeros(size(param.T));
l = zeros(size(param.T));
alpha_1 = zeros(size(param.T));

for k=1:param.T-1
    inputs = {x(:,k)};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    objective = solutions{3};
    if diagnostics == 1
        error('The problem is infeasible');
    end
    % Check if states are in tset
    mpc_in_tset(k) = Tset_mpt3.contains(x(:,k));
    lqr_in_tset(k) = Tset_mpt3.contains(x_lqr(:,k));    
    
    %update all the cost function to verify the assumption
    alpha_F = 1/2*x(:,k)'*cost.maxP*x(:,k);
    alpha_F_10 = 1/2*x(:,k)'*cost.maxP_10*x(:,k);
    obj(k)= objective;
    alpha_f(k) = alpha_F;
    alpha_f_10(k) = alpha_F_10;
    V_f(k) = 1/2*x(:,k)'*cost.Qf*x(:,k);
    l(k) =objective - 1/2*X(:,dim.N+1)'*cost.Qf*X(:,dim.N+1);
    alpha = 1/2*x(:,k)'*cost.minQ*x(:,k);
    alpha_1(k) = alpha;

    % Select the first input only
    u_rec(:,k) = U(1);
    u_lqr(:,k) = max(min(-LTI.K*x_lqr(:,k),con.umax),con.umin); % -LTI.K*x_lqr(:,k); 

    % Compute the state/output evolution
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u_rec(:,k); % + LTI.Bdist*d(:,k)
    x_lqr(:,k+1) = LTI.A*x_lqr(:,k) + LTI.B*u_lqr(:,k);

    % Check assumptions
    a(k) = Vf(x(:,k+1)) - Vf(x(:,k));
    b(k) = - ell(x(:,k), u_rec(:,k));
end
lqr_in_tset = ~lqr_in_tset; lqr_in_tset = ~lqr_in_tset;
mpc_in_tset = ~mpc_in_tset; mpc_in_tset = ~mpc_in_tset;

% plot_state(x, u_rec, 'State and Input Evolution with MPC linear state control');

figure(1);clf;
sgtitle(sprintf("LQR and MPC comparison (N = %d, equilibrium deviation = %.2f deg)", ...
    dim.N, 180/pi*param.eps))
subplot(3,1,1);

stairs(param.time, 180/(pi)*x_lqr(1,:), LineWidth=1.2); hold on;
stairs(param.time, 180/(pi)*x(1,:), LineWidth=1.2, LineStyle='-');
xlabel("Time (s)")
ylabel("x_1 (deg)")
grid on; hold off

subplot(3,1,2);
stairs(param.time, 180/(pi)*x_lqr(2,:), LineWidth=1.2); hold on;
stairs(param.time, 180/(pi)*x(2,:), LineWidth=1.2, LineStyle='-');
xlabel("Time (s)")
ylabel("x_2 (deg)")
grid on; hold off

subplot(3,1,3);
stairs(param.time, u_lqr(1,:), LineWidth=1.2); hold on;
stairs(param.time, u_rec(1,:), LineWidth=1.2, LineStyle='-'); 
xlabel("Time (s)")
ylabel("input (V)")
legend("LQR", "MPC", Location="best")
grid on; hold off

figure(2); clf;
stairs(param.time,a ,LineWidth=1.2); hold on;
stairs(param.time,b ,LineWidth=1.2, LineStyle="--");
xlabel("Time (s)")
legend({'$V_{f}(f(x,u)) - V_{f}(x)$', '$-\ell(x,u)$'}, ...
       'Interpreter', 'latex', ...
       'Location', 'best');

% Feasibility map for MPC
theta2_vals_deg    = -10:1:10;     % in degrees
theta2dot_vals_rad = -1.3:0.1:1.3;       % in rad/s

% Prepare a 3D matrix to store feasibility (true/false)
feasMap = false(length(theta2_vals_deg), length(theta2dot_vals_rad));%, length(theta2_vals_deg));

% Loop over all (theta2, theta2_dot, t
for i = 1:length(theta2_vals_deg)
   for j = 1:length(theta2dot_vals_rad)
            
       % Convert theta2 from deg to rad for internal calculations
       theta2_rad = theta2_vals_deg(i)*pi/180;
            
       % Build the initial state:
       % x0 = [theta1; theta2; theta1_dot; theta2_dot; i_motor]
       % with theta1=0, i_motor=0, and varying theta1_dot, theta2, theta2_dot
       x0 = [0; theta2_rad; 0; theta2dot_vals_rad(j); 0];
       inputs = {x0};
       [solutions,diagnostics] = controller{inputs};
            
       % If diagnostics == 0, then feasible
       feasMap(j, i) = (diagnostics == 0);
   end
end

% 3D Plot of the Feasibility Map
figure(3);clf; hold on; grid on;

for i = 1:length(theta2_vals_deg)
    for j = 1:length(theta2dot_vals_rad)

    % Get coordinate values
    xCoord = theta2_vals_deg(i)*pi/180;
    yCoord = theta2dot_vals_rad(j);
    % Pick color based on feasibility
    if feasMap(j, i)
        plot(xCoord, yCoord, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 3);
        else
           plot(xCoord, yCoord, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3);
        end
    end
end


view();hold on; % Set the 3D view
axis tight;


fixed_dims = [1,3,5];
fixed_vals = [0,0,0];

% Stati da fissare: x1, x3, xn = 0
S_mpt3 = Tset_mpt3.slice(fixed_dims,fixed_vals);

S_mpt3.plot('color', [0.678, 0.847, 0.902]);


xlabel('x_2 [rad]');
ylabel('x_4 [rad/s]');
sgtitle('Controllable state and Terminal Set for x_1 = x_3 = x_5 = 0');
grid on;

Tset_mpt3.isBounded()

%plot assumption 2.14(b), check initial condition in Xf or Xn
figure; clf;
stairs(param.time(1:end-1),log(V_f), 'LineWidth', 1.2); hold on;
stairs(param.time(1:end-1),log(alpha_f), 'LineWidth', 1.2);
stairs(param.time(1:end-1),log(alpha_1),'LineWidth', 1.2);
stairs(param.time(1:end-1),log(l), 'LineWidth', 2);
legend({'$V_f(x)$', '$\alpha_f(x)$', '$\alpha_1(x)$', '$\ell(x,u)$'}, ...
    'Interpreter', 'latex', 'Location', 'Best');

xlabel('Time(s)', 'Interpreter', 'latex');
ylabel('Log scale', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
box on;

%plot assumption 2.17 check initial condition in Xn
figure; clf;
stairs(param.time(1:end-1),log(obj), 'LineWidth', 1.2); hold on;
stairs(param.time(1:end-1),log(alpha_f_10), 'LineWidth', 1.2);
legend({'$V_N^{0}(x)$', '$\alpha(x)$'}, ...
    'Interpreter', 'latex', 'Location', 'Best');

xlabel('Time(s)', 'Interpreter', 'latex');
ylabel('Log scale', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
box on;

%% MPC Horizon tests
horizons = [5,30,100];
param.eps =4*pi/180;
param.time = 0:param.Ts:10;
param.T = length(param.time);    % simulation number of steps
% Initialize plot

figure(4); clf;
sgtitle(sprintf("Effect of Horizon on MPC (equilibrium deviation = %.2f deg)", ...
        180/pi*param.eps))

h = gcf;
for i = 1:length(horizons)
    N = horizons(i);
    % Initialize solver
    u = sdpvar(repmat(dim.nu,1,N),ones(1,N)); 
    x = sdpvar(repmat(dim.nx,1,N+1),ones(1,N+1));
    
    constraints = [];
    objective = 0;
    x_ref = [0;0;0;0;0];
    for k = 1:N
        objective = objective + 1/2*(x{k}-x_ref)'*cost.Q*(x{k}-x_ref) + 1/2*u{k}'*cost.R*u{k};
        constraints = [constraints, x{k+1} == LTI.A*x{k} + LTI.B*u{k}];
        constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmin <= x{k+1}<= con.xmax];
    end
    constraints = [constraints; con.Ff*x{N+1} <= con.ff];
    % constraints = [constraints; con.xmin <= x{N+1}<= con.xmax];
    
    objective = objective + 1/2*(x{N+1}-x_ref)'*cost.Qf*(x{N+1}-x_ref);
    
    parameters_in = x{1};
    solutions_out = {[u{:}], [x{:}]};
    
    controller = optimizer(constraints, objective,sdpsettings('solver','quadprog'),parameters_in,solutions_out);
    
    % 4. Simulation : Receding horizon implementation for the constrained control problem
    x = zeros(dim.nx,param.T);  % state vector with 0 at linearization point
    x(:,1) = [0;param.eps;0;0;0];  % starting from the up position of the pendulum with
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
    
        % Compute the state/output evolution
        x(:,k+1) = LTI.A*x(:,k) + LTI.B*u_rec(:,k); % + LTI.Bdist*d(:,k)
    end

    % plot results
    subplot(3,1,1); hold on;
    stairs(param.time, 180/(pi)*x(1,:), LineWidth=1.2);
    xlabel("Time (s)")
    ylabel("x_1 (deg)")
    grid on; hold off
    
    subplot(3,1,2); hold on;
    stairs(param.time, 180/(pi)*x(2,:), LineWidth=1.2);
    xlabel("Time (s)")
    ylabel("x_2 (deg)")
    grid on; hold off
    
    subplot(3,1,3); hold on;
    stairs(param.time, u_rec(1,:), LineWidth=1.2); 
    xlabel("Time (s)")
    ylabel("input (V)")
    grid on; hold off

end
subplot(3,1,3); hold on;
legend(sprintf("N = %d", horizons(1)),sprintf("N = %d", horizons(2)),sprintf("N = %d", horizons(3)))

%% MPC implementation - Non Linear dynamics; full state knowledge
clc;

param.eps = 4.968*pi/180;

disp("MPC implementation - Non Linear dynamics; full state knowledge")
u_nonlin = zeros(dim.nu,param.T); % input vector
u = u_nonlin;
x_nonlin = zeros(dim.nx, param.T);
x = x_nonlin;
x_lqr = x;
xref  = [0;0;0;0;0];

x_nonlin(:,1) = xref + [0;param.eps;0;0;0];
x(:,1) = [0;param.eps;0;0;0];
x_lqr(:,1) = x(:,1);

for k=1:(param.T-1)

    % Get the nonlinear control action
    inputs = {x_nonlin(:,k)-xref};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The non lin problem is infeasible');
    end
    u_nonlin(:,k) = U(1);
    
    % Get the linearized control action
    inputs = {x(:,k)};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The lin problem is infeasible');
    end
    u(:,k) = U(1);

    % sim real siyst
    % Simulate the nonlinear dynamics over one time step using ODE45
    [~, y_next] = ode45(@(t, y) furuta_nonlinear(y, u_nonlin(:,k), 0), [0 param.Ts], x_nonlin(:,k));
    x_nonlin(:,k+1) = y_next(end,:)';

    % Simulate the linear system
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u(:,k);

    ulqr =  max(min(-LTI.K*x_lqr(:,k),con.umax),con.umin);
    
    % Sim lqr on non lin
    [~, y_next] = ode45(@(t, y) furuta_nonlinear(y, ulqr , 0), [0 param.Ts], x_lqr(:,k));
    x_lqr(:,k+1) = y_next(end,:)';
end

figure(5);clf;
sgtitle(sprintf("Non linear and linearized dynamics MPC comparison \n (N = %d, equilibrium deviation = %.2f deg)", ...
    dim.N, 180/pi*param.eps))
subplot(3,1,1);

stairs(param.time, 180/(pi)*(x_nonlin(1,:)-xref(1)), LineWidth=1.2); hold on;
stairs(param.time, 180/(pi)*x(1,:), LineWidth=1.2, LineStyle='-');
% stairs(param.time, 180/(pi)*x_lqr(1,:), LineWidth=1.2, LineStyle='--');
xlabel("Time (s)")
ylabel("x_1 (deg)")
grid on; hold off

subplot(3,1,2);
stairs(param.time, 180/(pi)*(x_nonlin(2,:)-xref(2)), LineWidth=1.2); hold on;
stairs(param.time, 180/(pi)*x(2,:), LineWidth=1.2);
% stairs(param.time, 180/(pi)*x_lqr(2,:), LineWidth=1.2);

xlabel("Time (s)")
ylabel("x_2 (deg)")
grid on; hold off

subplot(3,1,3);
stairs(param.time, u_nonlin(1,:), LineWidth=1.2); hold on;
stairs(param.time, u(1,:), LineWidth=1.2); 

xlabel("Time (s)")
ylabel("input (V)")
legend("Non-Linear", "Linearized", Location="best")
grid on; hold off
