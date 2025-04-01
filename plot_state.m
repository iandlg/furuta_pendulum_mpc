function plot_state(x, u_rec, plot_title, x_hat)
    % PLOTTING RESULTS
    figure;
    clf;
    
    subplot(3,2,1);
    if nargin > 3 && ~isempty(x_hat)
        plot(x(1,:), '--', 'Color', 'black'); hold on;
        stairs(x_hat(1,:));
        legend('Real State', 'Estimated State');
    else
        stairs(x(1,:));
        legend('Real State');
    end
    title('State x_1 (\theta_1)');
    grid on;
    
    subplot(3,2,2);
    if nargin > 3 && ~isempty(x_hat)
        plot(x(2,:), '--', 'Color', 'black'); hold on;
        stairs(x_hat(2,:));
        legend('Real State', 'Estimated State');
    else
        stairs(x(2,:));
        legend('Real State');
    end
    title('State x_2 (\theta_2)');
    grid on;
    yline(0, '--r', 'Reference 0');  % Reference line for theta2
    
    subplot(3,2,3);
    if nargin > 3 && ~isempty(x_hat)
        plot(x(3,:), '--', 'Color', 'black'); hold on;
        stairs(x_hat(3,:));
        legend('Real State', 'Estimated State');
    else
        stairs(x(3,:));
        legend('Real State');
    end
    title('State x_3 (\theta_1 dot)');
    grid on;
    
    subplot(3,2,4);
    if nargin > 3 && ~isempty(x_hat)
        plot(x(4,:), '--', 'Color', 'black'); hold on;
        stairs(x_hat(4,:));
        legend('Real State', 'Estimated State');
    else
        stairs(x(4,:));
        legend('Real State');
    end
    title('State x_4 (\theta_2 dot)');
    grid on;
    
    subplot(3,2,5);
    if nargin > 3 && ~isempty(x_hat)
        plot(x(5,:), '--', 'Color', 'black'); hold on;
        stairs(x_hat(5,:));
        legend('Real State', 'Estimated State');
    else
        stairs(x(5,:));
        legend('Real State');
    end
    title('State x_5 (i)');
    grid on;
    
    subplot(3,2,6);
    stairs(u_rec(1,:));
    title('Input u (V)');
    legend('Control Input');
    grid on;
    
    xlabel('Time (s)');
    sgtitle(plot_title);
end
