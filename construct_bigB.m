function B_big = construct_bigB(A, B, N)
    % Get the dimensions of A and B
    [n, ~] = size(A);
    [~, m] = size(B);
    
    % Initialize B_big as a zero matrix
    B_big = zeros(n*(N+1), m*(N+1));
    
    % Construct the lower block triangular structure
    for i = 1:N+1
        for j = 1:i-1  % Ensuring that diagonal entries are zero
            B_big((i-1)*n+1:i*n, (j-1)*m+1:j*m) = A^(i-j-1) * B;
        end
    end
end
