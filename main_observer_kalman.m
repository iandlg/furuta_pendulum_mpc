% furuta Pendulum Simulation (Nonlinear Dynamics with Inputs)

clear; clc; close all;
load("mpt3_controller_data.mat");

%% MPC with observer, partial state knoledge

disp("MPC implementation - Observer and kalman")
% Get linearized system
param.y_eq = [0; 0; 0; 0; 0];  % Equilibrium point [theta1, theta2, theta1_dot, theta2_dot, i_motor]
param.u_eq = [0; 0];           % Equilibrium input [u1, u2]
param.Ts = 0.2;                %Sampling Time

[LTI.A, LTI.B, LTI.Bdist] = get_lin_dynamics(param.y_eq,param.u_eq,param.Ts);

% Simulation parameters
dim.N = 10;      % horizon
dim.nx = size(LTI.A,1);      % system order
dim.nu = 1;     % input order
dim.nd = 1;     % number of disturbance
param.time = 0:param.Ts:10;
param.T = length(param.time);    % simulation number of steps
cost.Q = diag([0.01;100;0.1;10;1]);  % Weighting on states (identity matrix)
cost.R = 1;             % Weighting on control input (scalar, as there's only one control input)
param.Anoise = 0.0001;
% Let Qf = P (solution to DARE)
[LTI.K, cost.Qf, ~] = dlqr(LTI.A, LTI.B, cost.Q, cost.R); % Optimal feedback gain
cost.Qf = 2*cost.Qf;

close all;clc;
disp("MPC implementation - Linear dynamics; Observer")
param.eps = 0.5*pi/180;
yref = [0;0];
x0 = [0;0;0;0;0];
LTI.x0 = x0;

dim.nd = 1;
LTI.C = [1 0 0 0 0;
         0 1 0 0 0]; % assume access to theta1 and 2 angles

dim.ny = size(LTI.C, 1);
LTI.Cdist = zeros(dim.ny,dim.nd);

LTIe.A = [LTI.A, LTI.Bdist; 
          zeros(dim.nd, dim.nx), eye(dim.nd)];

LTIe.B = [LTI.B; zeros(dim.nd,dim.nu)];

LTIe.C = [LTI.C,zeros(2,dim.nd)];


LTIe.x0 = [LTI.x0; 0];

LTIe.yref = yref;

dime.nx=dim.nx+dim.nd;     %state dimension
dime.nu=1;                 %input dimension
dime.ny=2;                 %output dimension
dime.N=10;                  %horizon

weighte.Q=blkdiag(cost.Q,2*ones(dim.nd));            %weight on output
weighte.R=cost.R;                                   %weight on input
weighte.P=blkdiag(cost.Qf,zeros(dim.nd));  

xe=zeros(dime.nx,param.T+1);
y=zeros(dime.ny,param.T+1);
u_rec_obs=zeros(dime.nu,param.T);
xehat=zeros(dime.nx,param.T+1);

xe(:,1)=LTIe.x0;
xehat(:,1)=[-0.1*pi/180; 0.1*pi/180; 0; 0; 0; 0];

measurement_noise = param.Anoise*normrnd(0,1,dim.ny,1);
y(:,1)=LTIe.C*LTIe.x0 +  measurement_noise;


 % L = place(LTIe.A',LTIe.C',[0.6; 0.55; 0.5;0.65 ;0.7; 0.4])';
 [L, ~, ~] = dlqr(LTIe.A', LTIe.C', weighte.Q, weighte.R); % Optimal feedback gain
 L = L';
options = sdpsettings('verbose',0,'solver','quadprog');

con.xmax = [100; 10*pi/180; 100; 100; 100];
con.xmin = -con.xmax;
con.xmaxe = [100; 10*pi/180; 100; 100; 100; 1];
con.xmine = -con.xmaxe;
con.umax = 10;
con.umin = -con.umax;
con.Ff = Tset_mpt3.A; con.ff = Tset_mpt3.b;


% optimizer for the controller
u = sdpvar(repmat(dime.nu,1,dime.N),ones(1,dime.N)); 
x = sdpvar(repmat(dime.nx,1,dime.N+1),ones(1,dime.N+1));
xr = sdpvar(dime.nx,1);
ur = sdpvar(dime.nu, 1);

constraints = [];
objective = 0;
for k = 1:dime.N
 objective = objective + 1/2*(x{k}-xr)'*weighte.Q*(x{k}-xr) + 1/2*(u{k}-ur)'*weighte.R*(u{k}-ur);
 constraints = [constraints, x{k+1} == LTIe.A*x{k} + LTIe.B*u{k}];
 constraints = [constraints, con.umin <= u{k}<= con.umax, con.xmine <= x{k+1}<= con.xmaxe];
end
constraints = [constraints; con.Ff*x{dim.N+1}(1:5) <= con.ff];

objective = objective + 1/2*(x{dime.N+1}-xr)'*weighte.P*(x{dime.N+1}-xr);

parameters_in = {x{1}, xr, ur};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,options,parameters_in,solutions_out);

for k=1:param.T
    if(k>=param.T/4&& k<=param.T/4+5) 
        xe(end,k) = 0.0008;
    else
        xe(end,k) = 0;
    end
    xe_0=xe(:,k);  
    dhat=xehat(end-dim.nd+1:end,k);
    
    [ineqconstraint, eqconstraints]= constraintsgen(LTI, dim, dhat, yref,con);
    [xref, uref] = optimalss(LTI, dim, cost, ineqconstraint, eqconstraints);
    xre = [xref; dhat];
    
     inputs = {xehat(:,k), xre, uref};
     [solutions,diagnostics] = controller{inputs};    
     U = solutions{1};
     X = solutions{2};
     if diagnostics == 1
         error('The problem is infeasible for controller');
     end     

    u_rec_obs(:,k)=U(1:dim.nu);
    % Compute the state/output evolution
    process_noise = param.Anoise * normrnd(0,1, [dime.nx, 1]);

    xe(:,k+1)=LTIe.A*xe_0 + LTIe.B*u_rec_obs(:,k) + process_noise;

    measurement_noise = param.Anoise*normrnd(0,1,dim.ny,1);
    y(:,k+1)=LTIe.C*xe(:,k+1)+ measurement_noise;
        
    % Update extended-state estimation

    xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec_obs(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));

end

% plot_state(xe , u_rec_obs , 'State Evolution with output MPC observer',xehat);

% Output MPC - with input and measurement noise, Kalman fiklter
disp("Output MPC - with input and measurement noise, Kalman filter")

param.eps = 1*pi/180; % deviation from equilibrium

Ob = obsv(LTI.A, LTI.C);
disp(['The rank of the observability matrix of the linearized matrix pair (A C) is : ', num2str(rank(Ob))])

% Kalman Filter variables
cov.pos = eye(dime.nx);
cov.measurement = eye(dime.ny);
cov.motion = 0.1*eye(dime.nx);


% Simulation : Receding horizon implementation for the constrained control problem
xhat = zeros(dime.nx,param.T);  % state vector with 0 at linearization point
xhat(:,1) = [-0.1*pi/180; 0.1*pi/180; 0; 0; 0; 0];  % starting from the up position of the pendulum with
x = zeros(dime.nx,param.T);
x(:,1) = [0;0;0;0;0;0];       % real starting position for state vector
u_rec_kal = zeros(dim.nu,param.T); % input vector

for k=1:param.T
    % Do measurement update 
    measurement_noise = param.Anoise*normrnd(0,1,dime.ny,1);
    y = LTIe.C*x(:,k) + measurement_noise;
    [xhat(:,k), cov.pos] = measurement_update(xhat(:,k), y, LTIe, cov);
    
    if(k>=param.T/4&& k<=param.T/4+5)
        x(end,k) = 0.0008;
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
    u_rec_kal(:,k) = U(1:dim.nu);

    % Compute the state/output evolution
    process_noise = param.Anoise*zeros(dime.nx, 1); process_noise(end,1) = param.Anoise* normrnd(0,1);
    x(:,k+1) = LTIe.A*x(:,k) + LTIe.B*u_rec_kal(:,k) + process_noise;
    
    % Do dynamic update for next k
    [xhat(:,k+1), cov.pos] = dynamic_update(xhat(:,k), u_rec_kal(:,k), LTIe, cov);
    
end
% Note we apply uniform noise to the measurement and noise only to the
% input


% plot_state(x, u_rec_kal, 'State Evolution with output MPC on linearized system,Kalman filter',xhat);

%
% Time vector
t = (0:param.T) * param.Ts;
figure;
% Theta1 subplot
subplot(3,1,1)
p1 = stairs(t, xe(1,:) * 180/pi, '-', 'Color', 'b'); hold on;
p2 = stairs(t, xehat(1,:) * 180/pi, '--', 'Color', 'b');
p3 = stairs(t, x(1,:) * 180/pi, '-', 'Color', 'r');
p4 = stairs(t, xhat(1,:) * 180/pi, '--', 'Color', 'r');
xlim([0,10]);
ylabel('$\theta_1$ [deg]', 'Interpreter', 'latex');
grid on;

legend([p1 p2 p3 p4], ...
    {'$\theta_1$ Real (Obs)', '$\theta_1$ Est. (Obs)', ...
     '$\theta_1$ Real (Kal)', '$\theta_1$ Est. (Kal)'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast');

% Theta2 subplot
subplot(3,1,2)
stairs(t, xe(2,:) * 180/pi, '-', 'Color', 'b'); hold on;
stairs(t, xehat(2,:) * 180/pi, '--', 'Color','b');
stairs(t, x(2,:) * 180/pi, '-', 'Color', 'r');
stairs(t, xhat(2,:) * 180/pi, '--', 'Color', 'r');
xlim([0,10]);
ylabel('$\theta_2$ [deg]', 'Interpreter', 'latex');
grid on;

% Input subplot
subplot(3,1,3)
p5 = stairs(t(1:end-1), u_rec_obs(1,:) , '-', 'Color', [0 0.6 0], 'LineWidth', 1.2); hold on;
p6 = stairs(t(1:end-1), u_rec_kal(1,:) , '-', 'Color', [0.6 0 0.6], 'LineWidth', 1.2);
xlim([0,10]);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Input [V]', 'Interpreter', 'latex');
grid on;

legend([p5 p6], {'Input (Observer)', 'Input (Kal)'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast');

% Global title
sgtitle('Observer vs. Kalman, Disturbance rejection without noise');