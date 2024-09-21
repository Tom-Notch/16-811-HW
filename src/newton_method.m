function root = newton_method(f, x0, tol, max_iter)
    % Initialize variables
    syms x real;

    f = str2sym(f);
    df = diff(f, x);

    f_handle = matlabFunction(f);
    df_handle = matlabFunction(df);

    disp("f_handle:")
    disp(f_handle);
    disp("df_handle:")
    disp(df_handle);

    x_current = double(x0);
    iter = 0;
    % Loop until the maximum number of iterations is reached
    while iter < max_iter
        % Calculate the derivative of the function at x
        % Calculate the next approximation of the root using Newton's method
        x_next = x_current - f_handle(x_current) / df_handle(x_current);
        % Check if the absolute difference between the current and next approximation is less than the tolerance
        if abs(f_handle(x_next)) < tol
            root = x_next;
            disp("Preliminary stop at iteration: " + iter);
            return;
        end

        % Update the current approximation and iteration count
        x_current = x_next;
        iter = iter + 1;
    end

    % If the maximum number of iterations is reached, return the current approximation
    root = x_current;
end
