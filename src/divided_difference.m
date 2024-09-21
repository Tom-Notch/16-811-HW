function polynomial_interpolation = divided_difference(X, F)
    % divided_difference constructs the interpolating polynomial
    % using Newton's Divided Difference method and returns the symbolic function.
    %
    % Inputs:
    %   X  - vector containing the x values: [x0, x1, ..., xn]
    %   F  - vector containing the corresponding function values: [f(x0), f(x1), ..., f(xn)]
    %
    % Output:
    %   polynomial_interpolation - the symbolic function representing the interpolating polynomial

    % Assert that X and F are of the same length
    assert(length(X) == length(F), 'X and F must have the same length.');

    syms x; % Define symbolic variable for x
    n = length(X); % Number of data points
    diff_table = zeros(n, n); % Initialize the divided difference table

    % Fill the first column with the function values
    diff_table(:, 1) = F(:); % Column vector of function values

    % Compute the divided differences
    for j = 2:n

        for i = 1:n - j + 1
            diff_table(i, j) = (diff_table(i + 1, j - 1) - diff_table(i, j - 1)) / (X(i + j - 1) - X(i));
        end

    end

    % Construct the interpolating polynomial symbolically
    polynomial_interpolation = diff_table(1, 1); % The first term (f(x0))
    product_term = 1; % To hold the product of (x - x0), (x - x1), ...

    for k = 2:n
        product_term = product_term * (x - X(k - 1)); % (x - x0)(x - x1)...(x - x(k-2))
        polynomial_interpolation = polynomial_interpolation + diff_table(1, k) * product_term;
    end

    % Simplify the polynomial
    polynomial_interpolation = simplify(polynomial_interpolation);
end
