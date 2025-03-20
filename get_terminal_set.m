function terminal_set = get_terminal_set(A,B,input_set, state_set,max_iter, tol)
% get_terminal_set If the function converges, the the resulting set is the maximal stabilizing set of the origin.
%                 It computes the terminal set
% Parameters : 
% Initialize :
[n,~] = size(A);
equilibrium_state = zeros(n,1);
terminal_set = Polyhedron('V',equilibrium_state');
disp(['Terminal set dimension at iteration ', num2str(0),' is ', num2str(terminal_set.Dim)]);


for i = 1:max_iter
    % Compute next terminal state
    B_U = input_set.affineMap(B);
    disp(['Affine map of U on B volume at iteration ', num2str(i),' is ', num2str(B_U.volume())]);
    disp(['Affine map of U on B dimension at iteration ', num2str(i),' is ', num2str(B_U.Dim)]);
    K_B_U = terminal_set.plus(-B_U);
    disp(['Minkovski K and B_U volume at iteration ', num2str(i),' is ', num2str(K_B_U.volume())]);
    K_B_U_A = K_B_U.invAffineMap(A);
    next_terminal_set = K_B_U_A.intersect(state_set);
    
    % Check for convergence
    set_difference = next_terminal_set.mldivide(terminal_set);
    volume = set_difference.volume();
    disp(['Converged at iteration ', num2str(i), ' with volume ', num2str(volume)]);
    % if volume < tol
    %     disp(['Converged at iteration ', num2str(i), ' with volume ', num2str(volume)]);
    %     break;
    % end
    terminal_set = next_terminal_set;
end
end
