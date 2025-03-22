function output = generate_reference(time,end_angle,should_plot)
% Define simulation time parameters
output = zeros(2, length(time));   % Initialize 2xN output matrix

% Define time segments
t1 = 3;                % First segment end (constant at 0)
t2 = 40;                % Second segment end (smooth transition to 2pi)
t3 = 50;                % Third segment end (constant at 2pi)

% Find indices corresponding to these times
idx1 = time <= t1;
idx2 = time > t1 & time <= t2;
idx3 = time > t2;

% First row: angle in radians
% First segment: constant at 0
output(1, idx1) = 0;

% Second segment: smooth transition from 0 to 2pi using sine acceleration profile
% This creates an S-curve for smooth acceleration and deceleration
t_norm = (time(idx2) - t1) / (t2 - t1);  % Normalized time (0 to 1)
% Using sine acceleration profile for smooth start and end
output(1, idx2) = end_angle * (t_norm - sin(2*pi*t_norm)/(2*pi));

% Third segment: constant at 2pi
output(1, idx3) = end_angle;

% Second row: all zeros
output(2, :) = pi;

if should_plot == true
    % Plot the results
    figure;
    subplot(2,1,1);
    plot(output(1,:), 'b-', 'LineWidth', 2);
    grid on;
    title('Angle (rad)');
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    
    subplot(2,1,2);
    plot(output(2,:), 'r-', 'LineWidth', 2);
    grid on;
    title('Second Parameter (always zero)');
    xlabel('Time (s)');
    ylabel('Value');
end
end
