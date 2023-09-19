function [outputArg1] = GenerateICs(num_functions,N,y_min,y_max)
%Generates smooth cubic spline initial data...
%N = number of spatial discretization points
%L = length of domain
%ymin = minimum interpolation value (function value can go slightly beyond)
%ymax = maximum interpolation value (function value can go slightly beyond)
Data = NaN(num_functions, N); %eventual output is num_func by N+1 array
%where i-th row corresponds to i-th initial condition.
% Generate random smooth functions for initial data
for i = 1:num_functions
    % Generate random control points
    control_points = unifrnd(y_min, y_max, 1, 3);

    % Add the first point as the last point to enforce periodicity
    control_points = [control_points, control_points(1)];

    % Fit a smooth curve through the control points using a cubic spline
    xvals = linspace(0, 1, 4);
    cs = csape(xvals, control_points, 'periodic');

    % Evaluate the curve on the given domain
    xvals = linspace(0, 1, N);
    y = fnval(cs, xvals);

    % Store the generated data in the Data array
    Data(i, :) = y;
end

outputArg1 = Data;

end