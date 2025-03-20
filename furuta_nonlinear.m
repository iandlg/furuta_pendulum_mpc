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
