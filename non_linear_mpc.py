import casadi as ca
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.linalg import cho_factor, cho_solve
from dataclasses import dataclass, field
import numpy as np
import scipy.io

def furuta_dynamics():
    # Dimensions
    dim_x = 5
    dim_u = 1

    # Symbolic variables
    q = ca.SX.sym('q', dim_x)
    u1 = ca.SX.sym('u', dim_u)

    # Extract state variables
    q1, q2, q3, q4, q5 = ca.vertsplit(q, 1)
    # u1, u2 = ca.vertsplit(u, 1)

    # Parameters
    L1 = 0.278  # Length parameter 1 [m]
    L2 = 0.300  # Length parameter 2 [m]
    l1 = 0.150  # Center of mass distance for first link [m]
    l2 = 0.148  # Center of mass distance for second link [m]
    m1 = 0.300  # Mass of first link [kg]
    m2 = 0.075  # Mass of second link [kg]
    J1 = 2.48e-2  # Inertia of first link [kg*m^2]
    J2 = 3.86e-3  # Inertia of second link [kg*m^2]
    b1 = 1.00e-4  # Damping coefficient for first link
    b2 = 2.80e-4  # Damping coefficient for second link
    g = 9.81  # Gravitational acceleration [m/s^2]
    Km = 0.090  # Motor constant
    Lm = 0.005  # Motor inductance [H]
    Rm = 7.80  # Motor resistance [Ohm]
    u2 = 0      # assume no disturbance for now

    J0_hat = J1 + m1 * l1**2 + m2 * L1**2
    J2_hat = J2 + m2 * l2**2

    # Dynamics
    q1_dot = q3

    q2_dot = q4

    q3_dot = (-J2_hat * b1 * q3 + m2 * L1 * l2 * ca.cos(q2) * b2 * q4 - J2_hat**2 * ca.sin(2 * q2) * q3 * q4 - 0.5 * J2_hat * m2 * L1 * l2 * ca.cos(q2) * ca.sin(2 * q2) * q3**2 + J2_hat * m2 * L1 * l2 * ca.sin(q2) * q4**2 + J2_hat * Km * q5 - m2 * L1 * l2 * ca.cos(q2) * u2 + 0.5 * m2**2 * l2**2 * L1 * ca.sin(2 * q2) * g) / (J0_hat * J2_hat + J2_hat**2 * ca.sin(q2)**2 - m2**2 * L1**2 * l2**2 * ca.cos(q2)**2)

    q4_dot = (m2 * L1 * l2 * ca.cos(q2) * b1 * q3 - b2 * (J0_hat + J2_hat * ca.sin(q2)**2) * q4 + m2 * L1 * l2 * J2_hat * ca.cos(q2) * ca.sin(2 * q2) * q3 * q4 - 0.5 * ca.sin(2 * q2) * (J0_hat * J2_hat + J2_hat**2 * ca.sin(q2)**2) * q3**2 - 0.5 * m2**2 * L1**2 * l2**2 * ca.sin(2 * q2) * q4**2 - m2 * L1 * l2 * ca.cos(q2) * Km * q5 + (J0_hat + J2_hat * ca.sin(q2)**2) * u2 - m2 * l2 * ca.sin(q2) * (J0_hat + J2_hat * ca.sin(q2)**2) * g) / (J0_hat * J2_hat + J2_hat**2 * ca.sin(q2)**2 - m2**2 * L1**2 * l2**2 * ca.cos(q2)**2)

    q5_dot = (u1 - Rm * q5 - Km * q3) / Lm

    # Concatenate derivatives
    q_dot = ca.vertcat(q1_dot, q2_dot, q3_dot, q4_dot, q5_dot)

    # Continuous dynamics function
    fc_furuta = ca.Function('fc_furuta', [q, u1], [q_dot])

    # RK4 integration
    dt = 0.05
    k1 = fc_furuta(q, u1)
    k2 = fc_furuta(q + dt/2 * k1, u1)
    k3 = fc_furuta(q + dt/2 * k2, u1)
    k4 = fc_furuta(q + dt * k3, u1)
    q_next = q + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # Discrete dynamics function
    f_rk4_furuta = ca.Function('f_rk4_furuta', [q, u1], [q_next])

    jac_dyn_q = ca.jacobian(q_next, q)
    jac_dyn_u = ca.jacobian(q_next, u1)
    jac_dyn_q_fun = ca.Function('jac_dyn_q_fun', [q, u1], [jac_dyn_q])
    jac_dyn_u_fun = ca.Function('jac_dyn_u_fun', [q, u1], [jac_dyn_u])

    return f_rk4_furuta, fc_furuta, jac_dyn_q_fun, jac_dyn_u_fun

def get_dynamics_matrices(x, u):
    f_rk4_furuta, fc_furuta, jac_dyn_q_fun, jac_dyn_u_fun = furuta_dynamics()
    A = jac_dyn_q_fun(x, u)
    B = jac_dyn_u_fun(x, u)
    return A.full(), B.full() # Convert to dense matrices

def solve_TO(x0, N, dt, x_ref, u_ref):
    # Get the dynamics
    f_rk4_furuta, _, _, _ = furuta_dynamics()
    
    # Dimensions
    dim_x = 5
    dim_u = 1

    # Cost weights
    Q = np.diag([1, 5000000, 1, 1, 1])
    R = np.array([[1]])
    Qf = 200 * Q
    
    # Define the state and control decision variables
    X = ca.SX.sym('x', dim_x, N + 1)
    U = ca.SX.sym('u', dim_u, N)
        
    # Quadratic cost function
    cost = 0.
    for k in range(N):
        x_err = X[:,k] - x_ref
        u_err = U[:,k] - u_ref 
        cost += ca.mtimes([x_err.T, Q, x_err]) + ca.mtimes([u_err.T, R, u_err])
    # Terminal cost
    x_err = X[:,-1] - x_ref
    cost += ca.mtimes([x_err.T, Qf, x_err])
    
    # Dynamics constraints
    g_dyn = [X[:,0] - x0]
    # RK4 integration
    for k in range(N):
        x_k = X[:,k]
        u_k = U[:,k]
        # RK4 step
        x_next = f_rk4_furuta(x_k, u_k)
        g_dyn.append(x_next - X[:,k+1])
    
    # Concatenate constraints
    g = ca.vertcat(*g_dyn)
    g_min = 0.
    g_max = 0.
    
    # Stack up the decision variables
    w = ca.vertcat(ca.reshape(X, -1, 1), ca.reshape(U, -1, 1))
    # Here, you can set bounds on states and controls if necessary
    # For simplicity, we assume no bounds (can be modified as needed)
    # Bounds on states and controls
    bound = ca.DM.ones(dim_x * (N + 1) + dim_u * N, 1)
    bound[:dim_x * (N + 1)+1, 0] = ca.inf
    bound[1+dim_x * (N + 1):, 0] = 10

    w_max = bound
    w_min = -bound
    # Create the NLP problem
    nlp = {'x': w, 'f': cost, 'g': g}
    # Create the solver instance
    opts = {
        'ipopt.print_level': 0,        # Disable IPOPT-specific output
        'print_time': False,           # Disable CasADi timing print
        'ipopt.sb': 'yes',             # Suppress IPOPT banner
        'verbose': False               # Disable any verbose output
    }
    solver = ca.nlpsol('solver', 'ipopt', nlp, 
                       opts # Specify the solver options 
    ) 
    
    X_init = np.zeros((dim_x, N + 1))
    U_init = np.zeros((dim_u, N))
    w_init = ca.vertcat(X_init.reshape((-1, 1)), U_init.reshape((-1, 1)))
    
    # Solve the NLP
    solution = solver(x0=w_init, lbx=w_min, ubx=w_max, lbg=g_min, ubg=g_max)
    
    # Extract the solution
    w_opt = solution['x']
    X_opt = w_opt[:dim_x * (N + 1)].reshape((dim_x, N + 1)).full()
    U_opt = w_opt[dim_x * (N + 1):].reshape((dim_u, N)).full()
    
    return X_opt, U_opt

def stage_cost(params, x, u, k):
    Q = params.Q
    R = params.R
    dx = x - params.xf
    du = u - params.uf  
    return 0.5 * dx.T @ Q @ dx + 0.5 * du.T @ R @ du

def final_cost(params, x):
    N = params.N
    Qf = params.Qf
    dx = x - params.xf
    return 0.5 * dx.T @ Qf @ dx

def stage_cost_expansion(params, x, u, k):
    Q = params.Q
    R = params.R
    dx = x - params.xf
    du = u - params.uf 
    # l_xx, l_ux, l_uu, l_x, l_u
    return deepcopy(Q), \
        np.zeros((R.shape[0], Q.shape[0])), \
        deepcopy(R), \
        Q @ dx, \
        R @ du

def final_cost_expansion(params, x):
    N = params.N
    Qf = params.Qf
    dx = x - params.xf
    # lf_xx, lf_x
    return deepcopy(Qf), \
        Qf @ dx

def trajectory_cost(params, x_trj, u_trj):
    N = params.N
    cost = 0.
    for k in range(N):
        cost += stage_cost(params, x_trj[k, :], u_trj[k, :], k)
    cost += final_cost(params, x_trj[N, :])
    return cost

def backward_pass(params, x_trj, u_trj, regu):
    symmetrize = lambda x: (x + x.T) / 2   
    N = params.N
    dim_x = x_trj.shape[1]
    dim_u = u_trj.shape[1]
        
    K_trj = np.zeros([N, dim_u, dim_x])
    d_trj = np.zeros([N, dim_u])
    expected_cost_redu = 0.

    # final cost expansion
    V_xx, V_x = final_cost_expansion(params, x_trj[N, :])

    for k in range(N-1, -1, -1):
        # dynamics jacobians
        A, B = get_dynamics_matrices(x_trj[k, :], u_trj[k, :])

        # stage cost expansion
        l_xx, l_ux, l_uu, l_x, l_u = stage_cost_expansion(params, x_trj[k, :], u_trj[k, :], k)

        # Q function expansion
        Q_x  = l_x + A.T @ V_x
        Q_u  = l_u + B.T @ V_x
        Q_xx = l_xx + A.T @ V_xx @ A
        Q_uu = l_uu + B.T @ V_xx @ B
        Q_ux = l_ux + B.T @ V_xx @ A

        # add regularization to ensure that Q_uu is invertible and well conditioned        
        Q_uu_regu = Q_uu + np.eye(dim_u) * regu
        Q_uu_regu = symmetrize(Q_uu_regu)

        chofact = cho_factor(Q_uu_regu)
        K = -cho_solve(chofact, Q_ux)
        d = -cho_solve(chofact, Q_u)
        K_trj[k, :, :] = K
        d_trj[k, :]    = d
        
        # cost-to-go
        V_xx = Q_xx + K.T @ Q_uu @ K + K.T @ Q_ux + Q_ux.T @ K
        V_xx = symmetrize(V_xx)
        V_x  = Q_x  + K.T @ Q_uu @ d + K.T @ Q_u  + Q_ux.T @ d

        # expected cost reduction
        expected_cost_redu += -Q_u.T @ d - 0.5 * d.T @ Q_uu @ d

    return K_trj, d_trj, expected_cost_redu

def forward_pass(params, 
                 x_trj, u_trj, 
                 K_trj, d_trj, 
                 cost,
                 f_dyn):
    N = params.N
    dim_x = x_trj.shape[1]
    dim_u = u_trj.shape[1]
    
    x_trj_new = np.zeros((N + 1, dim_x))
    x_trj_new[0, :] = x_trj[0, :]
    u_trj_new = np.zeros((N, dim_u))
    alpha = 1.
    
    # line search
    for _ in range(params.max_ls_iter):
        for k in range(N):
            u_trj_new[k, :] = u_trj[k, :] + K_trj[k, :, :] @ (x_trj_new[k, :] - x_trj[k, :]) + alpha * d_trj[k, :]
            x_trj_new[k+1, :] = f_dyn(x_trj_new[k, :], u_trj_new[k, :]).full().flatten()
        cost_new = trajectory_cost(params, x_trj_new, u_trj_new)

        if cost_new < cost:
            return x_trj_new, u_trj_new, alpha, cost_new
        alpha *= 0.5
    
    print('Line search failed!')
    alpha = 0.

    return x_trj, u_trj, alpha, cost

def update_regu(regu, regu_min, regu_max, alpha):
    if alpha == 0.:
        return min(regu_max, regu * 10)
    if alpha == 1.:
        return max(regu_min, regu * 0.1)
    return regu

def rollout(x0, u_trj, f_dyn):
    dim_x = x0.shape[0]
    N = u_trj.shape[0]
    x_trj = np.zeros((N + 1, dim_x))
    x_trj[0, :] = x0
    for k in range(N):
        x_trj[k+1, :] = f_dyn(x_trj[k, :], u_trj[k, :]).full().flatten()
    return x_trj

def run_ilqr(params, x0, f_dyn):
    N = params.N
    dim_u = params.dim_u
    u_trj = np.random.randn(N, dim_u) * 0.
    x_trj = rollout(x0, u_trj, f_dyn)
    # x_trj_hist = [x_trj]
    # u_trj_hist = [u_trj]
    
    regu = params.regu_init
    max_regu = params.max_regu
    min_regu = params.min_regu
    max_iter = params.max_iter
    
    for it in range(max_iter):
        traj_cost = trajectory_cost(params, x_trj, u_trj)
        # backward pass
        K_trj, d_trj, expected_cost_redu = backward_pass(params, x_trj, u_trj, regu)
        # forward pass
        x_trj_new, u_trj_new, alpha, traj_cost_new = forward_pass(params, x_trj, u_trj, K_trj, d_trj, traj_cost, f_dyn)
        x_trj = x_trj_new
        u_trj = u_trj_new
        regu = update_regu(regu, min_regu, max_regu, alpha)
        if alpha > 0 and np.abs(traj_cost_new - traj_cost) < params.cost_redu_tol:
            break
     
    return x_trj, u_trj

class OurParam:
    def __init__(self, N, dt):
        self.Q: np.ndarray = np.diag([1, 5000000, 1, 1, 1])
        self.R: np.ndarray = np.diag([1])
        self.Qf: np.ndarray = 200 * self.Q
        self.max_iter: int = 200
        self.regu_init: float = 0.001
        self.max_regu: float = 1e4
        self.min_regu: float = 0.001
        self.N : int = N
        self.dt : float = dt

        self.cost_redu_tol: float = 1e-5
        self.max_ls_iter: int = 20
        
        self.xf: np.ndarray = np.array([0., np.pi, 0., 0., 0.])
        self.uf: np.ndarray = np.array([0.])
        
        self.dim_x: int = 5
        self.dim_u: int = 1

        self.simu_iter: int = 400

        self.simu_time: np.ndarray = np.arange(0,self.simu_iter)*self.dt


if __name__ == "__main__":
    
    N = 250
    dt = 0.05
    x_ref = np.array([0, np.pi, 0, 0, 0])
    u_ref = np.array([0])
    x0 = np.array([0, 0, 0, 0, 0]) # np.pi+5*np.pi/180
    params = OurParam(N,dt)
    f_rk4_furuta, _, _, _ = furuta_dynamics()


    # Initialize variable
    X = np.zeros((params.dim_x, params.simu_iter),dtype=float)
    X[:,0] = x0
    U = np.zeros((params.dim_u, params.simu_iter),dtype=float)

    """    
    for k in range(params.simu_iter):
        print(k)
        x_curr = X[:,k]

        # Solve optimization
        _, U_opt_ca = solve_TO(x_curr, N, dt, x_ref, u_ref)

        # Get optimal control input
        u_opt = U_opt_ca[:,0]

        # Update the dynamics
        x_next = f_rk4_furuta(x_curr, u_opt)

        # Save the data for plotting
        X[:,k+1:k+2] = x_next
        U[:,k] = u_opt
    """
    # X_opt_ca, U_opt_ca = solve_TO(x0, N, dt, x_ref, u_ref)

    # X_opt_ilqr, U_opt_ilqr = run_ilqr(params, x0, f_rk4_furuta)

    X, U = solve_TO(x0, params.simu_iter, dt, x_ref, u_ref)

    time = np.arange(N + 1) * dt

    fig, axs = plt.subplots(3, 2, figsize=(10, 6))
    axs[0, 0].plot(params.simu_time, X[0,:-1])
    axs[0, 0].set(ylabel='q1 [m]')
    axs[0, 0].grid()

    axs[0, 1].plot(params.simu_time, X[1,:-1])
    axs[0, 1].set(ylabel='q2 [deg]')
    axs[0, 1].grid()

    axs[1, 0].plot(params.simu_time, X[2,:-1])
    axs[1, 0].set(xlabel='Time [s]', ylabel='q3 [m/s]')
    axs[1, 0].grid()

    axs[1, 1].plot(params.simu_time, X[3,:-1])
    axs[1, 1].set(xlabel='Time [s]', ylabel='q4 [rad/s]')
    axs[1, 1].grid()

    axs[2, 0].plot(params.simu_time, X[4,:-1])
    axs[2, 0].set(xlabel='Time [s]', ylabel='q5 [rad/s]')
    axs[2, 0].grid()    

    axs[2, 1].plot(params.simu_time, U[0,:])
    axs[2, 1].set(xlabel='Time [s]', ylabel='u [V]')
    axs[2, 1].grid()

    scipy.io.savemat("nonlin_mpc_results.mat",dict(X = X, U=U, N = N, Ts = dt, Q = params.Q, Qf = params.Qf, R = params.R)) 
    """
    fig, axs = plt.subplots(3, 2, figsize=(10, 6))
    axs[0, 0].plot(time, X_opt_ilqr[:,0])
    axs[0, 0].set(ylabel='q1 [m]')
    axs[0, 0].grid()

    axs[0, 1].plot(time, np.rad2deg(X_opt_ilqr[:,1]))
    axs[0, 1].set(ylabel='q2 [deg]')
    axs[0, 1].grid()

    axs[1, 0].plot(time, X_opt_ilqr[:,2])
    axs[1, 0].set(xlabel='Time [s]', ylabel='q3 [m/s]')
    axs[1, 0].grid()

    axs[1, 1].plot(time, X_opt_ilqr[:,3])
    axs[1, 1].set(xlabel='Time [s]', ylabel='q4 [rad/s]')
    axs[1, 1].grid()

    axs[2, 0].plot(time, X_opt_ilqr[:,4])
    axs[2, 0].set(xlabel='Time [s]', ylabel='q5 [rad/s]')
    axs[2, 0].grid()

    axs[2, 1].plot(time[:-1], U_opt_ilqr[:,0])
    axs[2, 1].set(xlabel='Time [s]', ylabel='u [rad/s]')
    axs[2, 1].grid()
    """
    

    plt.show()