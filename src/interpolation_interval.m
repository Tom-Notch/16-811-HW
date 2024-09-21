function min_interval = interpolation_interval(f, range, degree, error_tolerance)
    % interpolation_interval finds the minimum interval that satisfies the error tolerance to interpolate f with given degree of polynomial function.
    %
    % Inputs:
    %   f              - the function to interpolate
    %   range          - the range of x values: [a, b]
    %   degree         - the degree of the interpolating polynomial
    %   error_tolerance - the maximum error tolerance
    %
    % Output:
    %   min_interval - the minimum interval that satisfies the error tolerance

    syms x real; % Define symbolic variable for x

    f_sym = str2sym(f);

    % Derive coefficient function
    coefficient_function = diff(f_sym, x, degree + 1) / factorial(degree + 1);

    % Find the maximum value of the coefficient function in the range
    coefficient_function_handle = matlabFunction(-coefficient_function);
    max_coefficient_x = fminbnd(@(x) coefficient_function_handle(x), range(1), range(2));

    max_coefficient = subs(coefficient_function, max_coefficient_x);

    syms y h real;
    assumeAlso(h > 0 & h < 1);
    assumeAlso(y >= 0 & y <= degree * h);

    % Construct the error bound function
    error_bound_function = max_coefficient;

    for i = 0:degree
        error_bound_function = error_bound_function * (y - i * h);
    end

    error_bound_function = simplify(abs(error_bound_function));

    disp("error_bound_function:")
    disp(vpa(error_bound_function, 5));

    % Get the maximum expression of error_bound_function with y in range [0, degree * h] and h as positive constant
    % Can assume h will not be greater than 1
    % Resulting in a function with only h as variable
    derror_dy = diff(error_bound_function, y);

    disp("derror_dy:")
    disp(vpa(derror_dy, 5));

    critical_points = solve(derror_dy == 0, y, "Real", true);
    critical_points = critical_points(isAlways(critical_points >= 0 & critical_points <= degree * h));
    candidate_points = [0; critical_points; degree * h];

    disp("candidate_points:")
    disp(candidate_points);

    error_bound_function_values = subs(error_bound_function, y, candidate_points);

    % Find the maximum value and corresponding y
    max_error_bound_function = simplify(max(error_bound_function_values));

    disp("max_error_bound_function:")
    disp(vpa(max_error_bound_function, 5));

    % Solve for max_error_bound_function = error_tolerance / 2
    min_interval = solve(max_error_bound_function == (error_tolerance / 2), h, "Real", true);

end
