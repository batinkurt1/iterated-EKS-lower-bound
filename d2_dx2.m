function H = d2_dx2(f, x)
    % d2_dx2 computes the Hessian matrix of a scalar function f at point x
    %
    % Inputs:
    %   f - handle to the scalar-valued function, f: R^n -> R
    %   x - n-dimensional point at which to evaluate the Hessian
    %
    % Output:
    %   H - n x n Hessian matrix of second partial derivatives at x

    n = length(x);          % Number of variables
    H = zeros(n, n);       % Initialize Hessian matrix
    eps_step = 1e-3;       % Step size for finite differences

    for i = 1:n
        for j = i:n
            if i == j
                % Compute second derivative with respect to x_i twice
                x_forward = x;
                x_forward(i) = x_forward(i) + eps_step;
                
                x_backward = x;
                x_backward(i) = x_backward(i) - eps_step;
                
                f_forward = f(x_forward);
                f_backward = f(x_backward);
                f_current = f(x);
                
                H(i, i) = (f_forward - 2*f_current + f_backward) / (eps_step^2);
            else
                % Compute mixed partial derivatives ∂²f/∂x_i∂x_j
                x_pp = x;           % x + eps_step in both i and j
                x_pp(i) = x_pp(i) + eps_step;
                x_pp(j) = x_pp(j) + eps_step;
                
                x_p = x;            % x + eps_step in i
                x_p(i) = x_p(i) + eps_step;
                
                x_q = x;            % x + eps_step in j
                x_q(j) = x_q(j) + eps_step;
                
                f_pp = f(x_pp);
                f_p = f(x_p);
                f_q = f(x_q);
                f_current = f(x);
                
                % Mixed partial derivative using central difference
                H(i, j) = (f_pp - f_p - f_q + f_current) / (eps_step^2);
                H(j, i) = H(i, j); % Exploit symmetry
            end
        end
    end
end
