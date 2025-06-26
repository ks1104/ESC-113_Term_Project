function Y = rk4_solver(f, y0, x)
    % Custom Runge-Kutta 4th order solver
    % f: function handle for ODEs (dy/dt = f(t, y))
    % y0: initial values [S0; P0]
    % x: time vector

    n = length(x);    % Number of time steps
    h = x(2) - x(1);  % Step size (assuming uniform spacing)
    m = length(y0);   % Number of variables (e.g., S and P)

    Y = zeros(m, n);  % Preallocate solution matrix
    Y(:,1) = y0;      % Set initial conditions

    % Main RK4 loop
    for i = 1:n-1
        k1 = f(x(i), Y(:,i));                 % Slope at start
        k2 = f(x(i) + h/2, Y(:,i) + h*k1/2);  % Slope at midpoint 1
        k3 = f(x(i) + h/2, Y(:,i) + h*k2/2);  % Slope at midpoint 2
        k4 = f(x(i) + h,   Y(:,i) + h*k3);    % Slope at end

        % Weighted average of slopes
        Y(:,i+1) = Y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end
