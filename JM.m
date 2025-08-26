% init vars
maxTol = 1e-3;      
maxItr = 1e6; 
% set c & n values
cVal = [2,4]; 
nVal = [50, 100, 200, 400, 600, 800];

iterations_results = cell(1, length(cVal));

% for loop for Jacobi method
for c_id = 1:length(cVal)
    c = cVal(c_id);
    num_iterations = zeros(1, length(nVal));
    
    for i = 1:length(nVal)
        n = nVal(i);
        h = 1 / (n + 1);
        % make tridiag matrix A_n
        mDiag = (c / h^2) * ones(n, 1);
        offDiag = (-1 / h^2) * ones(n, 1);
        A_n = spdiags([offDiag mDiag offDiag], [-1 0 1], n, n);

        % r-side b_n
        b_n = ones(n, 1);

        % exac sol
        x_exact = A_n \ b_n;

        % init x0 
        x = zeros(n, 1);
        
        % get diag of A_n, then U and L too
        D = diag(diag(A_n));
        L = tril(A_n, -1);
        U = triu(A_n, 1);
        
        % get Dinv
        DInv = diag(1./diag(A_n));
        
        % get D^(-1) * b_n
        DInv_b = DInv * b_n;
        
        % get -D^(-1) * (L+U)
        B = -DInv * (L + U);
        
        % itr
        iter = 0;
        % error
        rel_error = 1;  
        
        while rel_error > maxTol && iter < maxItr
            iter = iter + 1;
            
            % comp new x using formula
            x_new = B * x + DInv_b;
          
            x = x_new;
            
            % comp rel error in infinity norm
            rel_error = norm(x - x_exact, inf) / norm(x_exact, inf);
        end
        
        % store # of itr
        num_iterations(i) = iter;
        
        % print the itr for each n
        % fprintf('for this c is %d and n is %d, then itr is %d\n', c, n, iter);
    end
    
    % store results
    iterations_results{c_id} = num_iterations;
    
    % plot for c[i]
    figure;
    plot(nVal, num_iterations, 'o-', 'LineWidth', 2);
    xlabel('Matrix Size (n)');
    ylabel('Number of Iterations');
    title(['Jacobi Iterations vs Matrix Size for c = ', num2str(c)]);
    grid on;
end

% log-log plot for big O
figure;
hold on;
colors = ['b', 'r'];
for c_id = 1:length(cVal)
    log_n = log(nVal);
    log_iters = log(iterations_results{c_id});
    coeffs = polyfit(log_n, log_iters, 1);
    % slope equal n in big O
    p = coeffs(1);  
    disp(['c = ', num2str(cVal(c_id)), ': O(n^', num2str(p), ')']);
    plot(log_n, log_iters, 'o-', 'Color', colors(mod(c_id - 1, length(colors)) + 1), ...
         'LineWidth', 2, 'DisplayName', ['c = ', num2str(cVal(c_id))]);
end
xlabel('log(n)');
ylabel('log(Iterations)');
title('Estimating Big-O Complexity (Jacobi Method)');
legend('Location', 'best');
grid on;
hold off;
